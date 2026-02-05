#!/usr/bin/env bash

# For all somatic loci in the "curated_somatic_sites" key of the config.yaml, confirm their presence in the somatic VCF.
# If absent, check for them in the pre-filtered VCF. If present in the pre-filtered dataset, add to the somatic VCF. 

# Load input arguments
threads=$1
somatic_vcf=$2
prefilter_vcf=$3
curated_sites=$4
n_sites=$5
tmp_out=$6
output=$7

# Test whether all sites are in somatic_vcf
n_somatic_sites=$(bcftools view -H -r ${curated_sites} ${somatic_vcf} | wc -l | xargs)
echo "Number of sites in somatic_vcf: ${n_somatic_sites}"

if [ "$n_sites" -eq "$n_somatic_sites" ]; then
    echo "All sites are in somatic_vcf"
    
    cp ${somatic_vcf} ${output}
    bcftools index \
        --threads ${threads} \
        --tbi ${output}

    # Create empty files as snakemake expects these temp(outputs)
    touch ${tmp_out}
    touch ${tmp_out}.tbi

    # Exit script without error
    exit 0

else
    echo "Not all sites are in somatic_vcf..."
    
    # Test whether all sites are in prefilter_vcf
    n_prefilter_sites=$(bcftools view -H -r ${curated_sites} ${prefilter_vcf} | wc -l | xargs)
    echo "Number of sites in prefilter vcf: ${n_prefilter_sites}"

    if [ "$n_prefilter_sites" -gt "$n_somatic_sites" ]; then
        echo "More sites found in prefilter_vcf"
        bcftools view \
            -r ${curated_sites} \
            -O z \
            -o ${tmp_out} \
            ${prefilter_vcf}

        bcftools index \
            --threads ${threads} \
            --tbi ${tmp_out}

        bcftools concat \
            --threads ${threads} \
            -a \
            -d all \
            -O z \
            -o ${output} \
            ${somatic_vcf} \
            ${tmp_out} 
        
        bcftools index \
            --threads ${threads} \
            --tbi ${output}
        
        exit 0

    else
        echo "No more sites found in prefilter_vcf"
        cp ${somatic_vcf} ${output}
        bcftools index \
            --threads ${threads} \
            --tbi ${output}

        # Create empty files as snakemake expects these temp(outputs)
        touch ${tmp_out}
        touch ${tmp_out}.tbi

        # Exit script without error
        exit 0
    fi
fi
