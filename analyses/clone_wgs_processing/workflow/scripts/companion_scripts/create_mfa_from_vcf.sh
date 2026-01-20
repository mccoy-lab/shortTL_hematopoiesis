#!/usr/bin/env bash

################################################################################################
# This script processes a gzipped VCF file containing biallelic SNPs for phylogenetic inference.
# It converts genotype codes to binary characters indicating the presence (HET or homALT), 
# absence (homREF), or missing status of the variant allele for each sample in the VCF.
#
# Usage: ./create_mfa_from_vcf.sh <input_vcf>
#
# Arguments:
#   input_vcf: A gzipped VCF file containing biallelic SNPs for phylogenetic inference.
#
# The script performs the following steps:
# 1. Checks if the input VCF file is provided and is gzipped.
# 2. For each sample in the VCF, it converts genotype codes to binary characters:
#    - '0' for homREF (reference homozygous)
#    - '1' for HET (heterozygous) or homALT (alternative homozygous)
#    - '-' for missing genotype
# 
# Output is printed to stdout
################################################################################################

input=$1 #vcf containing biallellic snps for phylogenetic inference

# check input length
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_vcf>"
    exit 1
fi

# check if input vcf is gzipped (crude suffix test) 
if [[ ${input} != *.gz ]]; then
    echo "Input VCF is not gzipped. Please provide a gzipped VCF file and associated index (tbi)."
    exit 1
fi

# for each sample in the vcf, convert genotype codes to binary characters of
# variant allele present (HET or homALT), absent (homREF), or missing. 
sample_list=`bcftools query --list-samples ${input}`
for sample in $sample_list; do
    echo ">${sample}"
    bcftools query --samples $sample -f "[%GT]," ${input} |\
        xargs |\
        tr -d '\n' |\
        sed 's/ //g' |\
        sed 's/,$/\n/' |\
        sed 's/|/\\/g' |\
        sed 's/0[\/|]0/0/g' |\
        sed 's/\.[\/|]\./-/g' |\
        sed 's/0[\/|]1/1/g' |\
        sed 's/1[\/|]0/1/g' |\
        sed 's/1[\/|]1/1/g' |\
        sed 's/,//g'
done