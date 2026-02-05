#!/usr/bin/env bash

# Remove samples that fail VAF binomial mixture model

#*Load inputs:
threads=$1
input_vcf=$2
sample_list=$3
output_tmp_vcf=$4
output_vcf=$5

#* Filter out bad VAF samples
echo "Filtering out samples failing VAF BMM..."
bcftools view \
    --threads ${threads} \
    -S ${sample_list} \
    --trim-alt-alleles \
    ${input_vcf} |\
bcftools view \
    --threads ${threads} \
    --include 'AC > 0' \
    -O z \
    -o ${output_tmp_vcf}

#* Index the filtered VCF
echo "Indexing filtered VCF..."
bcftools index \
    --threads ${threads} \
    --tbi ${output_tmp_vcf}

#* Recompute INFO field for final vcf output
echo "Recomputing AC, AN, AF, and DP INFO fields..."
bcftools annotate \
    --threads ${threads} \
    -x INFO/AC,INFO/AN,INFO/AF,INFO/DP \
    ${output_tmp_vcf} |\
    bcftools +fill-tags \
        --threads ${threads} -- \
            -t AN,AC,AF |\
    bcftools +fill-tags \
        --threads ${threads} \
        -O z \
        -o ${output_vcf} -- \
            -t 'DP:1=int(sum(FORMAT/DP))'

#* Index the final vcf output
echo "Indexing final VCF..."
bcftools index --threads ${threads} --tbi ${output_vcf}