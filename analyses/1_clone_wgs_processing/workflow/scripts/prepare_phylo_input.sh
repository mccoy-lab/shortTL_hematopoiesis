#!/usr/bin/env bash

# Sample germline reference for genotype calls and add to somatic VCF preior to phylogenetic reconstruction

# Load input
threads=$1
somatic_vcf=$2
germline_vcf=$3
germline_name=$4
tmp_germline_output=$5
output=$6

# Pull somatic variant loci from germline vcf
bcftools view \
    --samples ${germline_name} \
    -R ${somatic_vcf} \
    -O z \
    -o ${tmp_germline_output} \
    ${germline_vcf}
bcftools index \
    --tbi \
    ${tmp_germline_output}

# Merge germline vcf with somatic vcf and subset to snps
bcftools merge \
    -R ${somatic_vcf} \
    ${somatic_vcf} \
    ${tmp_germline_output} |\
bcftools view \
    --include 'TYPE="snp"' \
    -O v \
    -o ${output}