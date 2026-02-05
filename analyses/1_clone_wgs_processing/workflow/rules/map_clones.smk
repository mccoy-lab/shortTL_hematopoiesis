#!/usr/bin/env python3

"""
TITLE:
    Somatic variant calling for a population of blood clones.

AUTHORS:
    M. Tassia

DESCRIPTION:
    This pipeline is designed to call somatic variants for a population of blood clones.
    In addition to sequencing, mapping, variant calling stats, the major products of this
    pipeline are:
        1. A somatic VCF file with annotated variants
        2. A VCF of somatic SNPs for phylogenetic analysis

USAGE: 
    snakemake --configfile config/config.yaml --profile /path/to/profile_dir/ --use-conda -p # slurm
    snakemake --configfile config/config.yaml --use-conda --cores 16 -p # local

DETAILS:
    This pipeline assumes that all libraries are sequenced in a single run (i.e. same flowcell and run ID).
"""

#--- LIBRARY IMPORTS ---#
import pathlib
import pandas as pd

#--- CONFIGURATION ---#
# configure shell behavior for all rules
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail;")

#--- STATIC VARIABLE DECLARATION ---#

all_libraries = config['input']['libraries'] + ["consolidated_germline"]
chroms = list(range(1,23)) # define autosome list

#--- FUNCTION DECLARATION ---#
def get_sample_name(fp, id_str):
    """
    In a csv file (fp), use 'id_str' to match value in 
    SM_tag column and return the value in the column 'Subject_ID'.

    This function is used for renaming the raw fastq files after concatenating.
    """
    # Return "consolidated_germline" if id_str is "consolidated_germline"
    if id_str == "consolidated_germline":
        return "consolidated_germline"

    # read csv file into pandas dataframe
    df = pd.read_csv(fp) 

    # check if id_str has a suffix (e.g. _S1) and remove if present
    if "_" in id_str:
        id_str = id_str.split("_")[0]

    # return the value in the 'Subject_ID' column
    return str(df[df['SM_Tag'] == id_str]['Subject_ID'].values[0])

def get_sample_name_list(fp, id_list):
    """
    In a csv file (fp), search for each element in the input list (id_list),
    where each element in id_list will match to a value in the column SM_Tag,
    and return the values in Subject_ID as a list.
    """
    # read csv file into pandas dataframe
    df = pd.read_csv(fp)

    # create an empty list to store the sample names
    sample_names = []

    # iterate over each id in id_list
    for id_str in id_list:
        # check if id_str has a suffix (e.g. _S1) and remove if present
        if "_" in id_str:
            id_str = id_str.split("_")[0]
        
        # get the corresponding Subject_ID and append to the list
        sample_name = df[df['SM_Tag'] == id_str]['Subject_ID'].values[0]
        sample_names.append(sample_name)
    
    return sample_names

def get_sample_id(fp, name_str):
    """
    In a csv file (fp), use 'name_str' to match value in 
    Subject_ID column and return the value in the column 'SM_Tag'.

    This function is used for identifying germline samples.
    """
    # read csv file into pandas dataframe
    df = pd.read_csv(fp)

    # return the value in the 'SM_Tag' column
    return str(df[df['Subject_ID'] == name_str]['SM_Tag'].values[0])

def remove_low_coverage(fp, threshold):
    """
    Read a tsv file and return a list of samples with mean coverage above threshold.
    The tsv file is output from the generate_coverage_summaries rule.
    """
    # read tsv file into pandas dataframe
    df = pd.read_csv(fp, sep='\t')

    # filter samples with mean coverage below threshold
    return df[df['mean_cov'] >= threshold]['library'].tolist()

def convert_loci_to_list(curated_somatic_sites):
    """
    Read curated_somatic_sites as either a string or a list: 
        1. If string, return string.
        2. If list of length 1, return string.
        3. If list of length > 1, return comma-separated string.
    """
    if isinstance(curated_somatic_sites, str):
        return curated_somatic_sites
    elif isinstance(curated_somatic_sites, list):
        if len(curated_somatic_sites) == 1:
            return curated_somatic_sites[0]
        else:
            return ",".join(curated_somatic_sites)
    else:
        raise ValueError("curated_somatic_sites must be a string or a list of strings.")

#--- RULES ---#
localrules:
    aggregate_fastp_stats,
    plot_fastq_stats,
    aggregate_telseq_data,
    aggregate_bam_summary_stats,
    aggregate_bam_coverage_stats,
    generate_coverage_summaries,
    aggregate_mosdepth_global_dist,
    aggregate_mosdepth_summaries,
    genotyping_sensitivity_stats,
    somatic_filter,
    somatic_vcf_stats

rule fastp_stats:
# compute fastp stats for each consolidated library
    input:
        R1=config['input']['directory']+"/{library_lane}_R1_001.fastq.gz",
        R2=config['input']['directory']+"/{library_lane}_R2_001.fastq.gz"
    output: 
        html = config['output_dir'] + "/stats/fastp/{library_lane}_fastp.html",
        json = config['output_dir'] + "/stats/fastp/{library_lane}_fastp.json"
    conda:
        "../envs/stats_env.yaml"
    threads: 
        config['run_params']['fastp']['threads']
    resources:
        time = "4:00:00",
        mem_mb = config['run_params']['fastp']['mem_mb']
    shell:
        """
        fastp -V -i {input.R1} -I {input.R2} --json {output.json} --html {output.html} --report_title {wildcards.library_lane} --thread {threads}
        """

rule aggregate_fastp_stats:
# compile summary reports for all fastp results
    input: 
        json = expand(config['output_dir'] + "/stats/fastp/{sample}_{lane}_fastp.json",
                            sample = config['input']['libraries'],
                            lane = ["L001", "L002", "L003", "L004", "L005", "L006", "L007", "L008"]),
        aggregate_fastp = "workflow/scripts/read_json_stats.py"
    log: 
        stderr = "logs/aggregate_fastp_stats.log"
    output:
        tsv = config['output_dir'] + "/stats/aggregated_fastp_stats.tsv"
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        python {input.aggregate_fastp} {input.json} > {output.tsv} 2> {log.stderr}
        """

rule plot_fastq_stats:
# generate a plot that summarizes the per-sample read counts across runs and lanes
    input:
        tsv = rules.aggregate_fastp_stats.output.tsv,
        plot_script = "workflow/scripts/plot_fastq_stats.R"
    output:
        jpg = config['output_dir'] + "/stats/sample_read_counts.jpg"
    log:
        stderr = "logs/plot_fastq_stats.log"
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        Rscript {input.plot_script} {input.tsv} {output.jpg} 2> {log.stderr}
        """

rule map_libraries:
# map reads to reference genome per lane using bwa (with read grouping)
    input:
        R1=config['input']['directory']+"/{sample}_{lane}_R1_001.fastq.gz",
        R2=config['input']['directory']+"/{sample}_{lane}_R2_001.fastq.gz",
        ref = config['references']['grch38'],
        bwa_script = "workflow/scripts/bwa.sh"
    output:
        bam = temp(config['output_dir'] + "/bams/{sample}_{lane}.bam")
    wildcard_constraints:
        sample = "|".join(config['input']['libraries'])
    threads:
        config['run_params']['bwa']['threads']
    resources:
        mem_mb = config['run_params']['bwa']['mem_mb'],
        time = "24:00:00"
    conda:
        "../envs/analysis_env.yaml"
    shell:
        """
        bash {input.bwa_script} \
            {input.R1} \
            {input.R2} \
            {input.ref} \
            {output.bam} \
            {threads} \
            {wildcards.sample} \
            {wildcards.sample} 
        """

rule process_bams:
# sort, index, mark duplicates 
    input:
        bams = expand(config['output_dir'] + "/bams/{{sample}}_{lane}.bam",
                        lane = ["L001", "L002", "L003", "L004", "L005", "L006", "L007", "L008"])
    output: 
        processed_bam = config['output_dir'] + "/bams/{sample}.markdup.sort.bam"
    wildcard_constraints:
        sample = "|".join(config['input']['libraries'])
    params:
        out_dir=config['output_dir']
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['bwa']['mem_mb'],
        time = "48:00:00"
    conda:
        "../envs/analysis_env.yaml"
    shell:
        """
        #* Create temporary directory for intermediate files
        mkdir -p {params.out_dir}/tmp 

        #* Create input string for merging BAMs
        for file in {input.bams}; do
            bam_inputs+="-I $file "
        done

        #* Merge BAMs
        gatk MergeSamFiles \
            --java-options "-Xmx{resources.mem_mb}m" \
            -SO queryname \
            --USE_THREADING true \
            $bam_inputs \
            -O {params.out_dir}/tmp/{wildcards.sample}.sorted.bam

        #* Mark duplicates (assuming sequenced on Illumina NovaSeq)
        gatk MarkDuplicatesSpark \
            -I {params.out_dir}/tmp/{wildcards.sample}.sorted.bam \
            -O {output.processed_bam} \
            --optical-duplicate-pixel-distance 2500
        
        #* Index the processed BAM
        samtools index \
            -@ {threads} \
            {output.processed_bam}
        
        #* Remove intermediate files
        rm {params.out_dir}/tmp/{wildcards.sample}.sorted.bam
        """

rule telseq:
# run telseq to estimate somatic telomere length per sample
    input:
        bam = rules.process_bams.output.processed_bam,
    output:
        telseq = temp(config['output_dir'] + "/telseq/{sample}.telseq.txt")
    params:
        read_len = 150,
        repeat_thresh = 12,
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['telseq']['threads']
    resources:
        mem_mb = config['run_params']['telseq']['mem_mb'],
        time = "12:00:00"
    shell:
        """
        telseq -r {params.read_len} \
            -k {params.repeat_thresh} \
            -o {output.telseq} \
            {input.bam}
        """

rule aggregate_telseq_data:
# aggregate all telseq output (per sample) into a single tsv
   input:
       telseq_data = expand(config['output_dir'] + "/telseq/{sample}.telseq.txt",
                           sample = config['input']['libraries']),
   output:
       tsv = config['output_dir'] + "/stats/aggregated_telseq_data.tsv"
   conda:
       "../envs/stats_env.yaml"
   shell:
       """
       head -n 1 {input.telseq_data[0]} > {output.tsv}
       cat {input.telseq_data} |\
          grep -v 'LENGTH_ESTIMATE' >> {output.tsv}
       """

rule extract_pileup_window:
# extract pileup data for window specified in config['input']['pileup_window']
   input:
       bam = rules.process_bams.output.processed_bam,
       ref = config['references']['grch38']
   output: 
       pileup = temp(config['output_dir'] + "/pileups/{sample}.pileup.bam"),
       index = temp(config['output_dir'] + "/pileups/{sample}.pileup.bam.bai")
   params:
       region = config['input']['pileup_window']
   conda:
       "../envs/analysis_env.yaml"
   threads:
       config['run_params']['samtools']['threads']
   resources:
       mem_mb = config['run_params']['samtools']['mem_mb'],
       time = "1:00:00"
   shell:
       """
       #* Subset bams to focal region
       samtools view \
          -@ {threads} \
          -h \
          -b \
          -T {input.ref} \
          -o {output.pileup} \
          {input.bam} \
          {params.region}

        #* Index the pileup BAM
       samtools index \
          -@ {threads} \
          -b \
          {output.pileup}
       """

rule merge_pileup_windows:
# merge pileup data across samples
   input:
       pileups = expand(config['output_dir'] + "/pileups/{sample}.pileup.bam",
                       sample = config['input']['libraries']),
       pileup_indices = expand(config['output_dir'] + "/pileups/{sample}.pileup.bam.bai",
                       sample = config['input']['libraries'])
   output:
       merged_pileup = config['output_dir'] + "/bams/focal_pileup_window.merged.bam",
       merged_index = config['output_dir'] + "/bams/focal_pileup_window.merged.bam.bai"
   conda:
       "../envs/analysis_env.yaml"
   threads:
       config['run_params']['samtools']['threads']
   resources:
       mem_mb = config['run_params']['samtools']['mem_mb'],
       time = "1:00:00"
   shell:
       """
       #* Merge bams 
       samtools merge \
           -@ {threads} \
           -r \
           {output.merged_pileup} \
           {input.pileups}
        
        #* Index bams
        samtools index \
           -@ {threads} \
           -b \
           {output.merged_pileup}
        """

rule consolidated_germline_mapping:
# map germline sample reads to reference genome per germline sample+lane using bwa (with read grouping)
    input:
        R1=config['input']['directory']+"/{sample}_{lane}_R1_001.fastq.gz",
        R2=config['input']['directory']+"/{sample}_{lane}_R2_001.fastq.gz",
        ref = config['references']['grch38'],
        bwa_script = "workflow/scripts/bwa.sh"
    output:
        bam = temp(config['output_dir'] + "/bams/{sample}_{lane}.conslidated_germline.bam")
    wildcard_constraints:
        sample = "|".join(config['input']['germline_libraries'])
    threads:
        config['run_params']['bwa']['threads']
    resources:
        mem_mb = config['run_params']['bwa']['mem_mb'],
        time = "24:00:00"
    conda:
        "../envs/analysis_env.yaml"
    shell:
        """
        bash {input.bwa_script} \
            {input.R1} \
            {input.R2} \
            {input.ref} \
            {output.bam} \
            {threads} \
            consolidated_germline \
            consolidated_germline 
        """

rule consolidated_germline_processing:
# merge, sort, index, mark duplicates for consolidated_germline
    input:
        bams = expand(config['output_dir'] + "/bams/{sample}_{lane}.conslidated_germline.bam",
                        sample = config['input']['germline_libraries'],
                        lane = ["L001", "L002", "L003", "L004", "L005", "L006", "L007", "L008"])
    output: 
        processed_bam = config['output_dir'] + "/bams/consolidated_germline.markdup.sort.bam"
    wildcard_constraints:
        sample = "|".join(config['input']['germline_libraries'])
    params:
        sample_name="consolidated_germline",
        out_dir=config['output_dir']
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['bwa']['mem_mb'],
        time = "48:00:00"
    conda:
        "../envs/analysis_env.yaml"
    shell:
        """
        #* Create temporary directory for intermediate files
        mkdir -p {params.out_dir}/tmp 

        #* Create input string for merging BAMs
        for file in {input.bams}; do
            bam_inputs+="-I $file "
        done

        #* Merge BAMs
        gatk MergeSamFiles \
            --java-options "-Xmx{resources.mem_mb}m" \
            -SO queryname \
            --USE_THREADING true \
            $bam_inputs \
            -O {params.out_dir}/tmp/{params.sample_name}.sorted.bam
        
        #* Mark duplicates (assuming sequenced on Illumina NovaSeq)
        gatk MarkDuplicatesSpark \
            -I {params.out_dir}/tmp/{params.sample_name}.sorted.bam \
            -O {output.processed_bam} \
            --optical-duplicate-pixel-distance 2500

        #* Index the processed BAM
        samtools index \
            -@ {threads} \
            {output.processed_bam}

        #* Remove intermediate files
        rm {params.out_dir}/tmp/{params.sample_name}.sorted.bam
        """

rule bam_stats:
# generate bam stats
    input:
        bam = rules.process_bams.output.processed_bam,
        ref = config['references']['grch38']
    output:
        stats = config['output_dir'] + "/stats/bam/{sample}.bam.stats"
    threads:
        config['run_params']['samtools']['threads']
    resources:
        mem_mb = config['run_params']['bwa']['mem_mb'],
        time = "6:00:00"
    conda:  
        "../envs/analysis_env.yaml"
    shell:
        """
        samtools stats -@ {threads} \
            --reference {input.ref} \
            {input.bam} > {output.stats}
        """

rule aggregate_bam_summary_stats:
# compile summary reports for all bam stats
    input: 
        bam_stats = expand(config['output_dir'] + "/stats/bam/{sample}.bam.stats",
                            sample = all_libraries),
        aggregate_bam = "workflow/scripts/aggregate_bam_summary_stats.py"
    output:
        tsv = config['output_dir'] + "/stats/aggregated_bam_summary_stats.tsv"
    log:
        stderr = "logs/aggregate_bam_summary_stats.log"
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        python {input.aggregate_bam} {input.bam_stats} > {output.tsv} 2> {log.stderr}
        """

rule aggregate_bam_coverage_stats:
# compile coverage data for all bam stats
    input: 
        bam_stats = expand(config['output_dir'] + "/stats/bam/{sample}.bam.stats",
                            sample = all_libraries),
        aggregate_coverage = "workflow/scripts/aggregate_bam_coverage_stats.py"
    output:
        tsv = config['output_dir'] + "/stats/aggregated_bam_coverage_stats.tsv"
    log:
        stderr = "logs/aggregate_bam_coverage_stats.log"
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        python {input.aggregate_coverage} {input.bam_stats} > {output.tsv} 2> {log.stderr}
        """

rule generate_coverage_summaries:
# generate a plot and tsv that summarized the per-sample mean coverage
    input:
        coverage = rules.aggregate_bam_coverage_stats.output.tsv,
        summary_script = "workflow/scripts/summarize_coverage.R"
    output:
        tsv = config['output_dir'] + "/stats/coverage_summary.tsv",
        plot = config['output_dir'] + "/stats/coverage_summary.jpg"
    log:
        stderr = "logs/coverage_summary.log"
    params:
        ind = config['input']['individual']
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        Rscript {input.summary_script} {input.coverage} {params.ind} {output.tsv} {output.plot} 2> {log.stderr}
        """

rule mosdepth:
# run mosdepth to estimate chromosome-level coverage distributions
    input:
        bam = rules.process_bams.output.processed_bam
    output:
        summary = temp(config['output_dir'] + "/stats/mosdepth/{sample}.mosdepth.summary.txt"),
        depth_dist = temp(config['output_dir'] + "/stats/mosdepth/{sample}.mosdepth.global.dist.txt")
    params:
        out_dir = config['output_dir'] + "/stats/mosdepth"
    conda:
        "../envs/stats_env.yaml"
    threads:
        config['run_params']['mosdepth']['threads']
    resources:
        mem_mb = config['run_params']['mosdepth']['mem_mb'],
        time = "2:00:00"
    shell:
        """
        mosdepth -t {threads} -x -n {params.out_dir}/{wildcards.sample} {input.bam}
        """

rule aggregate_mosdepth_summaries:
    input:
        summaries = expand(config['output_dir'] + "/stats/mosdepth/{sample}.mosdepth.summary.txt",
                            sample = all_libraries),
        summary_script = "workflow/scripts/aggregate_mosdepth_summaries.py"
    output:
        tsv = config['output_dir'] + "/stats/mosdepth_summary.tsv"
    log:
        stderr = "logs/aggregate_mosdepth_summaries.log"
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        python {input.summary_script} {input.summaries} > {output.tsv} 2> {log.stderr}
        """

rule aggregate_mosdepth_global_dist:
    input:
        depth_dist = expand(config['output_dir'] + "/stats/mosdepth/{sample}.mosdepth.global.dist.txt",
                            sample = all_libraries),
        summary_script = "workflow/scripts/aggregate_mosdepth_global_dist.py",
        plot_script = "workflow/scripts/plot_mosdepth_depth_dist.R"
    output:
        tsv = config['output_dir'] + "/stats/mosdepth_depth_dist.tsv",
        jpg = config['output_dir'] + "/stats/mosdepth_depth_dist.jpg"
    log:
        stderr = "logs/aggregate_mosdepth_global_dist.log"
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        python {input.summary_script} {input.depth_dist} > {output.tsv} 2> {log.stderr}
        Rscript {input.plot_script} {output.tsv} {output.jpg} 2>> {log.stderr}
        """

rule create_reference_index:
# create fasta indices for haplotype caller if they do not exist
    input:
        ref = config['references']['grch38']
    output:
        index = config['references']['grch38'] + ".fai",
        dictionary = config['references']['grch38'].replace(".fa", ".dict")
    conda:
        "../envs/analysis_env.yaml"
    resources:
        mem_mb = config['run_params']['gatk']['mem_mb'],
        time = "2:00:00"
    shell:
        """
        samtools faidx {input.ref}
        gatk CreateSequenceDictionary \
            --java-options "-Xmx{resources.mem_mb}m" \
            -R {input.ref} -O {output.dictionary}
        """

rule haplotypecaller:
# HaplotypeCaller per chromosome per sample
    input:
        bam = rules.process_bams.output.processed_bam,
        ref = config['references']['grch38'],
        ref_index = rules.create_reference_index.output.index,
        ref_dict = rules.create_reference_index.output.dictionary
    output:
        gVCF = temp(config['output_dir'] + "/gVCFs/{sample}.chr{chrom}.g.vcf.gz"),
        index = temp(config['output_dir'] + "/gVCFs/{sample}.chr{chrom}.g.vcf.gz.tbi")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['gatk']['mem_mb'],
        time = "24:00:00"
    shell:
        """
        gatk HaplotypeCaller --java-options "-Xmx{resources.mem_mb}m" \
            --native-pair-hmm-threads {threads} \
            --emit-ref-confidence GVCF \
            --reference {input.ref} \
            --input {input.bam} \
            --intervals chr{wildcards.chrom} \
            --output {output.gVCF}
        """

rule genomicDBimport:
# Make genomeDB across sample gVCFs
    input:
        gvcfs = expand(config['output_dir'] + "/gVCFs/{sample}.chr{{chrom}}.g.vcf.gz",
                        sample = all_libraries),
        indices = expand(config['output_dir'] + "/gVCFs/{sample}.chr{{chrom}}.g.vcf.gz.tbi",
                        sample = all_libraries)
    output:
        db = temp(directory(config['output_dir'] + "/tmp/" + config['input']['individual'] + "_GenDB.chr{chrom}"))
    params:
        tmp_dir = config['output_dir'] + "/tmp",
        namemap = temp(config['output_dir'] + "/tmp/" + config['input']['individual'] + "_GenDB.chr{chrom}.namemap.txt")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['gatk']['mem_mb'],
        time = "48:00:00"
    shell:
        """
        #* Delete namemap if it exists
        if [ -f {params.namemap} ]; then
            rm {params.namemap}
        fi

        #* Create empty namemap for appending
        touch {params.namemap}

        #* Create namemap for GenomicsDBImport
        for file in {input.gvcfs}; do
            name=$(bcftools query --list-samples $file | xargs)
            printf "$name\t$file\n" >> {params.namemap}
        done

        #* Import gVCFs into GenomicsDB using namemape
        gatk GenomicsDBImport \
            --java-options "-Xmx{resources.mem_mb}m" \
            --sample-name-map {params.namemap} \
            --genomicsdb-workspace-path {output.db} \
            -L chr{wildcards.chrom} \
            --tmp-dir {params.tmp_dir} \
            --reader-threads {threads}

        #* Remove namemap
        rm {params.namemap}
        """

rule genotype_gvcfs:
# Using the genomicsDB, genotype samples
    input:
        db = rules.genomicDBimport.output.db,
        ref = config['references']['grch38']
    output:
        vcf = temp(config['output_dir'] + "/vcfs/" + config['input']['individual'] + "_chr{chrom}.vcf.gz"),
        index = temp(config['output_dir'] + "/vcfs/" + config['input']['individual'] + "_chr{chrom}.vcf.gz.tbi")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['gatk']['mem_mb'],
        time = "48:00:00"
    shell:
        """
        gatk GenotypeGVCFs \
            --java-options "-Xmx{resources.mem_mb}m" \
            -R {input.ref} \
            -V gendb://{input.db} \
            -O {output.vcf}
        """
    
rule vqsr:
# run VQSR per chromosome
    input:
        vcf = rules.genotype_gvcfs.output.vcf,
        index = rules.genotype_gvcfs.output.index,
        axiom = config['references']['axiom'],
        dbsnp = config['references']['dbsnp'],
        mills = config['references']['mills'],
        hapmap = config['references']['hapmap'],
        omni = config['references']['omni'],
        tgp = config['references']['tgp'],
        script = "workflow/scripts/VQSR.sh"
    output:
        vcf = temp(config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + "_chr{chrom}.vqsr.vcf.gz"),
        index = temp(config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + "_chr{chrom}.vqsr.vcf.gz.tbi")
    params:
        temp_dir = config['output_dir'] + "/tmp"
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['gatk']['mem_mb'],
        time = "24:00:00"
    shell:
        """
        mkdir -p {params.temp_dir}

        bash {input.script} -i {input.vcf} \
            -v {output.vcf} \
            -a {input.axiom} \
            -d {input.dbsnp} \
            -m {input.mills} \
            -h {input.hapmap} \
            -o {input.omni} \
            -t {input.tgp}  \
            -p {params.temp_dir} \
            -c {wildcards.chrom} \
            -x {resources.mem_mb}
        """

rule combine_VQSR_vcfs:
# combine each chromosomal vcf after VQSR
    input:
        vcfs = expand(config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + "_chr{chrom}.vqsr.vcf.gz",
                        chrom = chroms),
        indices = expand(config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + "_chr{chrom}.vqsr.vcf.gz.tbi",
                        chrom = chroms)
    output:
        vcf = config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + ".vqsr.vcf.gz",
        index = config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + ".vqsr.vcf.gz.tbi"
    params:
        tmp_dir = config['output_dir'] + "/tmp",
        input_list = temp(config['output_dir'] + "/tmp/input.list")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['gatk']['mem_mb'],
        time = "24:00:00"
    shell:
        """
        #* Create temporary directory for intermediate files
        mkdir -p {params.tmp_dir}

        #* Delete input list if it exists
        if [ -f {params.input_list} ]; then
            rm {params.input_list}
        fi

        #* Create input list for appending in loop below
        touch {params.input_list}

        #* Append each chromosomal VCF to the input list
        for file in {input.vcfs}; do
            printf "$file\n" >> {params.input_list}
        done

        #* Merge VCFs
        gatk MergeVcfs \
            --java-options "-Xmx{resources.mem_mb}m" \
            -I {params.input_list} \
            -O {output.vcf}

        #* Remove input list
        rm {params.input_list}
        """

rule fix_gatk_death_of_the_dot:
# recode genotypes as missing if DP=0 or SMPL_MAX(PL) = 0 (i.e. no reads support any genotype)
    input: 
        vcf = rules.combine_VQSR_vcfs.output.vcf
    output:
        vcf = config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + ".vqsr.recodedDOT.vcf.gz",
        index = config['output_dir'] + "/vqsr_vcf/" + config['input']['individual'] + ".vqsr.recodedDOT.vcf.gz.tbi"
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "6:00:00"
    shell:
        """
        bcftools +setGT --threads {threads} {input.vcf} -- -t q -i 'FORMAT/DP=0' -n . |\
            bcftools +setGT --threads {threads} -- -t q -i 'SMPL_MAX(PL) = 0' -n . |\
            bcftools annotate --threads {threads} -x INFO/AC,INFO/AN,INFO/AF,INFO/DP |\
            bcftools +fill-tags --threads {threads} -- -t AN,AC,AF |\
            bcftools +fill-tags --threads {threads} -O z -o {output.vcf} -- -t 'DP:1=int(sum(FORMAT/DP))'
        bcftools index --threads {threads} --tbi {output.vcf}
        """

rule genotyping_stats:
# compute variant calling stats for each sample using `bcftools stats` after removing VQSR data
    input:
        vcf = rules.fix_gatk_death_of_the_dot.output.vcf
    output:
        stats = config['output_dir'] + "/stats/genotyping_stats.txt"
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "1:00:00"
    shell:
        """
        bcftools stats -f PASS \
            --threads {threads} \
            --samples - \
            {input.vcf} > {output.stats}
        """

rule remove_positive_controls:
# remove sequencing positive control sample(s) from vcf if present
    input: 
        vcf = rules.fix_gatk_death_of_the_dot.output.vcf
    output:
        vcf = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.clones.vcf.gz"),
        index = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.clones.vcf.gz.tbi"),
        sample_list = temp(config['output_dir'] + "/filtered_vcfs/pos_controls.list")
    conda:
        "../envs/analysis_env.yaml"
    params:
        pos_controls_names = config['input']['germline_libraries'] + ["consolidated_germline"],
        pos_controls = config['input']['germline_libraries'] + ["consolidated_germline"]
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "6:00:00"
    shell:
        """
        #* QC printing to stdout
        echo {params.pos_controls}
        echo {params.pos_controls_names}

        #* Write positive control sample names to a file
        for sample in {params.pos_controls_names}; do
            echo "$sample" >> {output.sample_list}
        done

        #* Remove positive control samples from VCF and apply PASS filter
        bcftools view \
            --threads {threads} \
            -S ^{output.sample_list} \
            --trim-alt-alleles \
            -f PASS \
            {input.vcf} |\
        bcftools view \
            --threads {threads} \
            --include 'AC > 0' \
            -O z \
            -o {output.vcf} 
        
        #* Index the filtered VCF
        bcftools index \
            --threads {threads} \
            --tbi \
            {output.vcf}
        """

rule isolate_consolidated_germline_vcf:
# isolate the germline control (matched normal) from the parent VCF
    input: 
        vcf = rules.fix_gatk_death_of_the_dot.output.vcf
    output:
        vcf = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.consolidated_germline.vcf.gz",
        index = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.consolidated_germline.vcf.gz.tbi",
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "6:00:00"
    shell:
        """
        #* Extract consolidated germline sample from VCF and apply PASS filter
        bcftools view \
            --threads {threads} \
            -s consolidated_germline \
            -f PASS \
            -O z \
            -o {output.vcf} \
            {input.vcf}

        #* Index the filtered VCF
        bcftools index \
            --threads {threads} \
            --tbi \
            {output.vcf}
        """

rule prepare_sample_sensitivity:
# compute sensitivity per sample; output a tsv and a temporary vcf of high-confidence germline sites
    input:
        germ = rules.isolate_consolidated_germline_vcf.output.vcf,
        clones = rules.remove_positive_controls.output.vcf
    output:
        stats = temp(config['output_dir'] + "/stats/" + config['input']['individual'] + "_germline_HET_sites.tsv"),
        bcftools_stats = temp(config['output_dir'] + "/tmp/bcftools_stats.txt"),
        germline_het = temp(config['output_dir'] + "/tmp/germline_het.vcf.gz"),
        germline_het_index = temp(config['output_dir'] + "/tmp/germline_het.vcf.gz.tbi"),
        clone_germ = temp(config['output_dir'] + "/tmp/clone_germ.vcf.gz"),
        clone_germ_index = temp(config['output_dir'] + "/tmp/clone_germ.vcf.gz.tbi")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "2:00:00"
    shell:
        """
        #* Filter consolidated_germline variants to only high-confidence + biallelic sites
        echo "Preparing germline HET sites..."
        bcftools view \
            --threads {threads} \
            -m 2 \
            -M 2 \
            -a \
            {input.germ} |\
        bcftools view \
            -f PASS \
            --include 'COUNT(FMT/GT=="het" & FMT/GQ > 20) > 0' \
            -O z \
            -o {output.germline_het}
        bcftools index \
            --threads {threads} \
            --tbi {output.germline_het}
        
        #* Subset clonal data to only the high-confidence germline sites
        echo "Pulling germline HET sites from clones..."
        bcftools merge \
            --threads {threads} \
            -R {output.germline_het} \
            -O z \
            -o {output.clone_germ} \
            {input.clones} {output.germline_het}
        bcftools index \
            --threads {threads} \
            --tbi {output.clone_germ}

        #* Compute stats and subset to only per-sample stats (PSC block)
        echo "Computing stats on germline HET sites..."
        bcftools stats \
            --threads {threads} \
            --samples - \
            {output.clone_germ} > {output.bcftools_stats}
        grep 'PSC' {output.bcftools_stats} |\
            grep -v 'Per-sample counts' |\
            sed "s/^# PSC/PSC/" > {output.stats}
        """

rule genotyping_sensitivity_stats:
# generate a plot that summarizes the genotyping stats and sensitivity
    input:
        stats = rules.genotyping_stats.output.stats,
        sensitivity = rules.prepare_sample_sensitivity.output.stats,
        plot_script = "workflow/scripts/plot_clone_genotyping_stats.R"
    output:
        sensitivity = config['output_dir'] + "/stats/clone_sensitivity.tsv",
        jpg = config['output_dir'] + "/stats/clone_sensitivity.jpg",
        genotyping_psc = temp(config['output_dir'] + "/tmp/genotyping_psc.tsv")
    log:
        stderr = "logs/genotyping_sensitivity_stats.log"
    conda:
        "../envs/stats_env.yaml"
    shell:
        """
        #* Subset to only per-sample stats (PSC block)
        grep 'PSC' {input.stats} |\
            grep -v 'Per-sample counts' |\
            sed "s/^# PSC/PSC/" > {output.genotyping_psc}
        
        #* Generate clonal sensitivity tsv and plot
        Rscript {input.plot_script} \
            {output.genotyping_psc} \
            {input.sensitivity} \
            {output.sensitivity} \
            {output.jpg} 2> {log.stderr}
        """

rule filter_samples_by_sensitivity:
# Omit samples under sensitivty threshold defined in `config['run_params']['sensitivity_threshold']`
    input: 
        sensitivity = rules.genotyping_sensitivity_stats.output.sensitivity,
        vcf = rules.remove_positive_controls.output.vcf,
        index = rules.remove_positive_controls.output.index,
        filter_script = "workflow/scripts/filter_tsv_by_values.py"
    output:
        vcf = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.hs_clones.vcf.gz"),
        index = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.hs_clones.vcf.gz.tbi"),
        sample_list = temp(config['output_dir'] + "/filtered_vcfs/high_sensitivity_samples.list")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "3:00:00"
    shell:
        """
        #* Write high sensitivity samples to a file
        python {input.filter_script} \
            --tsv_path {input.sensitivity} \
            --output_path {output.sample_list} \
            --output_field sample \
            --filter_field sensitivity \
            --filter_threshold 0.5 \
            --mode greater_equal \
            --remove_consolidated_germline
            
        #* Filter out low sensitivity samples and remove sites without any observed HET genotypes
        bcftools view \
            --threads {threads} \
            -S {output.sample_list} \
            --trim-alt-alleles \
            -f PASS \
            {input.vcf} |\
        bcftools view \
            --threads {threads} \
            --include 'AC > 0' \
            -O z \
            -o {output.vcf} 

        #* Index the filtered VCF
        bcftools index \
            --threads {threads} \
            --tbi {output.vcf}
        """

rule recompute_vcf_fields:
# recompute AN,AC,AF fields, add VAF and DIST fields after removing samples
    input:
        vcf = rules.filter_samples_by_sensitivity.output.vcf,
        index = rules.filter_samples_by_sensitivity.output.index
    output:
        vcf = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.hs_clones.recomputedINFO.vcf.gz",
        index = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".recodedDOT.hs_clones.recomputedINFO.vcf.gz.tbi"
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "3:00:00"
    shell:
        """
        bcftools annotate --threads {threads} {input.vcf} -x INFO/AC,INFO/AN,INFO/AF,INFO/DP |\
        bcftools +fill-tags --threads {threads} -- -t AN,AC,AF,VAF |\
        bcftools +fill-tags --threads {threads} -- -t 'DP:1=int(sum(FORMAT/DP))' |\
        bcftools +variant-distance -O z -o {output.vcf}
        
        bcftools index --threads {threads} --tbi {output.vcf}
        """

rule somatic_filter_prep:
# prepare VCF for somatic filtering by removing sites:
#   1. Present in GnomAD (≥1% AF),
#   2. Present in 1KGP genomic accessibility mask,
#   3. Present in matched normal high-confidence sites,
#   4. Multiallelic,
#   5. Within 10bp of another variant
    input:
        vcf = rules.recompute_vcf_fields.output.vcf,
        index = rules.recompute_vcf_fields.output.index,
        germline_vcf = rules.prepare_sample_sensitivity.output.germline_het,
        germline_index = rules.prepare_sample_sensitivity.output.germline_het_index,
        gnomad = config['references']['gnomad'],
        mask = config['references']['tgp_accessibility'],
        script = "workflow/scripts/somatic_filter_prep.sh"
    output:
        vcf = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic_filter_prep.vcf.gz"),
        index = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic_filter_prep.vcf.gz.tbi")
    params:
        dist = config['run_params']['dist_threshold']
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "6:00:00"
    shell:
        """
        bash {input.script} \
            {input.vcf} \
            {input.germline_vcf} \
            {input.mask} \
            {input.gnomad} \
            {params.dist} \
            {threads} \
            {output.vcf}
        """

rule somatic_filter:
# Apply somatic filtering to the vcf
    input: 
        vcf = rules.somatic_filter_prep.output.vcf,
        index = rules.somatic_filter_prep.output.index,
        script = "workflow/scripts/SomaticFilterAmalgam.R"
    output:
        vcf = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic_filtered.vcf.gz"),
        index = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic_filtered.vcf.gz.tbi"),
        tmp = temp(config['output_dir'] + "/tmp/somatic_filter.vcf")
    log:
        "logs/somatic_filter.log"
    conda:
        "../envs/somatic_filter_env.yaml"
    threads:
        config['run_params']['somatic_filter']['threads']
    shell:
        """
        Rscript {input.script} \
            -i {input.vcf} \
            -o {output.tmp} \
            -c {threads} &> {log}
        bcftools view -O z \
            -o {output.vcf} \
            --threads {threads} \
            {output.tmp} &> {log}
        bcftools index --threads {threads} --tbi {output.vcf} &> {log}
        """

rule curate_somatic_sites:
# For variants known to be somatic, confirm their presence in the filtered vcf
# If absent, check pre-filtered vcf to see if the site was removed.
    input:
        somatic_vcf = rules.somatic_filter.output.vcf,
        somatic_index = rules.somatic_filter.output.index,
        prefilter_vcf = rules.recompute_vcf_fields.output.vcf,
        prefilter_index = rules.recompute_vcf_fields.output.index,
        script = "workflow/scripts/curate_somatic_sites.sh"
    output:
        vcf = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic.curated.vcf.gz"),
        index = temp(config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic.curated.vcf.gz.tbi"),
        tmp_out = temp(config['output_dir'] + "/tmp/curated_somatic_sites.vcf.gz"),
        tmp_index = temp(config['output_dir'] + "/tmp/curated_somatic_sites.vcf.gz.tbi")
    params:
        curated_sites = convert_loci_to_list(config['input']['curated_somatic_sites']),
        n_sites = len(config['input']['curated_somatic_sites'])
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "1:00:00"
    shell:
        """
        echo "Curated sites: {params.curated_sites}"
        echo "N. curated of sites: {params.n_sites}"

        bash {input.script} \
            {threads} \
            {input.somatic_vcf} \
            {input.prefilter_vcf} \
            {params.curated_sites} \
            {params.n_sites} \
            {output.tmp_out} \
            {output.vcf}
        """


rule somatic_vcf_stats:
# test evidence of clonal mixture,
# create file containing somatic vcf stats, 
# and plot somatic vcf stats
    input:
        vcf = rules.curate_somatic_sites.output.vcf,
        index = rules.curate_somatic_sites.output.index,
        sensitivity = rules.genotyping_sensitivity_stats.output.sensitivity,
        script = "workflow/scripts/plot_somatic_vaf_stats.R"
    output:
        jpg = config['output_dir'] + "/stats/somatic_vcf_stats.jpg",
        stats = config['output_dir'] + "/stats/somatic_vcf_stats.txt",
        bmm = config['output_dir'] + "/stats/somatic_VAF_bmm.tsv",
        psc = temp(config['output_dir'] + "/tmp/somatic_psc.tsv"),
        query_stats = temp(config['output_dir'] + "/tmp/somatic_stats.tsv"),
        allele_frequency = temp(config['output_dir'] + "/tmp/somatic_af.tsv")
    conda:
        "../envs/stats_env.yaml"
    log:
        "logs/somatic_vcf_stats.log"
    threads:
        config['run_params']['bcftools']['threads']
    shell:
        """
        #* Extract somatic allele frequencies
        echo "Extracting somatic allele frequencies..." &> {log}
        printf "site\tAF\n" > {output.allele_frequency}
        bcftools query \
           -f '%CHROM:%POS\t%AF\n' \
            {input.vcf} >> {output.allele_frequency} 2>> {log}

        #* Extract VAF and DP for each somatic variant per sample
        echo "Extracting VAF and DP per variant/sample..." &>> {log}
        printf "site\tsample\tDP\tVAF\n" > {output.query_stats}
        bcftools query \
            -f '[%CHROM:%POS\t%SAMPLE\t%DP\t%VAF\n]' \
            {input.vcf} >> {output.query_stats} 2>> {log}

        #* Run bcftools stats on curated somatic vcf
        echo "Running bcftools stats on somatic VCF..." &>> {log}
        bcftools stats \
            --samples - \
            {input.vcf} > {output.stats} 2>> {log}

        #* Extract PSC from bcftools stats output for plotting
        echo "Extracting PSC from bcftools stats output..." &>> {log}
        grep 'PSC' {output.stats} |\
            grep -v 'Per-sample counts' |\
            sed "s/^# PSC/PSC/" > {output.psc} 2>> {log}

        #* Run BMM on somatic VAF and plot curated somatic VCF stats
        echo "Running BMM on somatic VAFs..." &>> {log}
        Rscript {input.script} \
            {output.allele_frequency} \
            {output.query_stats} \
            {input.sensitivity} \
            {output.psc} \
            {output.bmm} \
            {output.jpg} 
        """

rule filter_by_vaf_bmm:
# filter by VAF BMM and recompute INFO stats
    input: 
        bmm_tsv = rules.somatic_vcf_stats.output.bmm,
        vcf = rules.curate_somatic_sites.output.vcf,
        index = rules.curate_somatic_sites.output.index,
        vaf_filter = "workflow/scripts/filter_by_vaf.sh",
        tsv_filter = "workflow/scripts/filter_tsv_by_values.py"
    output: 
        vcf = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic.vcf.gz",
        index = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic.vcf.gz.tbi",
        tmp_vcf = temp(config['output_dir'] + "/tmp/tmp.vcf.gz"),
        tmp_vcf_index = temp(config['output_dir'] + "/tmp/tmp.vcf.gz.tbi"),
        sample_list = temp(config['output_dir'] + "/filtered_vcfs/good_vaf_samples.list")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    resources:
        mem_mb = config['run_params']['bcftools']['mem_mb'],
        time = "3:00:00"
    shell:
        """
        #* Write VAF BMM-passing samples to sample list
        echo "Generating sample list..."
        python {input.tsv_filter} \
            --tsv_path {input.bmm_tsv} \
            --output_path {output.sample_list} \
            --output_field sample \
            --filter_field p \
            --filter_threshold 0.40 \
            --mode greater \
            --remove_consolidated_germline

        #* Filter out bad VAF samples and update INFO fields
        echo "Subsetting VCF and recoding INFO fields..."
        bash {input.vaf_filter} \
             {threads} \
             {input.vcf} \
             {output.sample_list} \
             {output.tmp_vcf} \
             {output.vcf}
        """

rule annotate_vcf:
# annotate somatic variants with funcotator
    input:
        vcf = rules.filter_by_vaf_bmm.output.vcf,
        index = rules.filter_by_vaf_bmm.output.index,
        ref = config['references']['grch38'],
        funcotator_resource = config['references']['funcotator_resource']
    output:
        vcf = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic.annotated.vcf.gz",
        index = config['output_dir'] + "/filtered_vcfs/" + config['input']['individual'] + ".somatic.annotated.vcf.gz.tbi",
        funcotator_resource = temp(directory(config['output_dir'] + "/tmp/funcotator"))
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['gatk']['threads']
    resources:
        mem_mb = config['run_params']['gatk']['mem_mb'],
        time = "6:00:00"
    shell:
        """
        #* Create funcotator resource directory for unpacking tarball
        mkdir -p {output.funcotator_resource}

        #* Extract input.funcotator_resource to output.funcotator_resource
        tar -xzvf {input.funcotator_resource} \
              -C {output.funcotator_resource} \
              --strip-components=1
              
        #* Annotate somatic variants with funcotator
        gatk Funcotator \
            --java-options "-Xmx{resources.mem_mb}m" \
            --ref-version hg38 \
            --data-sources-path {output.funcotator_resource} \
            --output-file-format VCF \
            --variant {input.vcf} \
            --reference {input.ref} \
            --output {output.vcf}
        """

rule sigprofiler:
# run SigProfilerAssigment to assign COSMIC somatic mutation signatures across samples
    input:
        vcf = rules.filter_by_vaf_bmm.output.vcf,
        index = rules.filter_by_vaf_bmm.output.index,
        script = "workflow/scripts/sigprofiler.py"
    output:
        sigprofiler_dir = directory(config['output_dir'] + "/sigprofiler"),
        metadata = config['output_dir'] + "/sigprofiler/JOB_METADATA_SPA.txt",
        vcf_split_dir = temp(directory(config['output_dir'] + "/sigprofiler/vcf_split")),
    params:
        genome = "GRCh38"
    conda:
        "../envs/sigprofiler_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    shell:
        """
        #* Create directory for sigprofiler input
        mkdir -p {output.vcf_split_dir}

        #* Split VCF into individual sample VCFs
        for sample in `bcftools query --list-samples {input.vcf}`; do
            bcftools view \
                --threads {threads} \
                -s $sample \
                {input.vcf} |\
            bcftools view \
                --threads {threads} \
                --include 'AC > 0' \
                -O v \
                -o {output.vcf_split_dir}/$sample.vcf
        done

        #* Run SigProfilerAssignment
        python {input.script} \
          {output.vcf_split_dir} \
          {output.sigprofiler_dir} \
          {params.genome}
        """

rule prep_phylogeny:
# merge germline genotypes for somatic variant calls for CellPhy input
    input:
        somatic_vcf = rules.filter_by_vaf_bmm.output.vcf,
        somatic_index = rules.filter_by_vaf_bmm.output.index,
        germline_vcf = rules.fix_gatk_death_of_the_dot.output.vcf,
        germline_index = rules.fix_gatk_death_of_the_dot.output.index,
        script = "workflow/scripts/prepare_phylo_input.sh"
    output:
        vcf = config['output_dir'] + "/phylogeny/" + config['input']['individual'] + ".snps_for_phylogeny.vcf",
        tmp_vcf = temp(config['output_dir'] + "/tmp/phylo_input.vcf.gz"),
        tmp_index = temp(config['output_dir'] + "/tmp/phylo_input.vcf.gz.tbi")
    conda:
        "../envs/analysis_env.yaml"
    threads:
        config['run_params']['bcftools']['threads']
    shell:
        """
        bash {input.script} \
            {threads} \
            {input.somatic_vcf} \
            {input.germline_vcf} \
            consolidated_germline \
            {output.tmp_vcf} \
            {output.vcf}
        """

        