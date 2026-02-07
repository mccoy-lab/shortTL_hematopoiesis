This directory contains scripts used for post-processing steps after completion of the 
full snakemake pipeline. The contents are described below:

script: cellphy.slurm.sh
description: Run CellPhy on the output from `rule prep_phylogeny`
details: copy this script to the directory containing the output from `rule prep_phylogeny`
    and replace <input.vcf> with the <name>.snps_for_phylogeny.vcf file contained in
    that directory.
usage: sbatch cellphy.slurm.sh

script: create_mfa_from_vcf.sh
description: Convert <name>.snps_for_phylogeny.vcf into a multiple fasta alignment of 
    binary states.
details: copy this script to the directory containing the output from `rule prep_phylogeny`.
    Run `./create_mfa_from_vcf.sh <input_vcf>`, where <input_vcf> is the file
     <name>.snps_for_phylogeny.vcf contained in that directory.
usage: ./create_mfa_from_vcf.sh <input_vcf>

script: fasta_to_phylip.py
description: Convert a multifasta alignment to phylip format
details:
usage: python fasta_to_phylip.py <input_fasta> <output_phylip>
