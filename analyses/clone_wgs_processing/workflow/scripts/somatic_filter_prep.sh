#!/usr/bin/env bash

# Run fixed filters on VCFs in preparation for the algorithmic filters

#Load inputs:
in_vcf=$1
in_germline=$2
in_mask=$3
in_gnomad=$4
dist=$5
threads=$6
out_vcf=$7

bcftools view \
	--threads ${threads} \
	-m 2 \
	-M 2 \
	-f PASS \
	${in_vcf} |\
bcftools filter \
	--threads ${threads} \
	--mask-file ${in_mask} \
	-s UnreliableSite -m + |\
bcftools filter \
	--threads ${threads} \
	--mask-file ${in_gnomad} \
	-s CommonPopVariant -m + |\
bcftools filter \
	--threads ${threads} \
	--mask-file ${in_germline} \
	-s GermlineHet -m + |\
bcftools view \
	-f PASS \
	--threads ${threads} |\
bcftools view \
	--threads ${threads} \
	--exclude "INFO/DIST <= ${dist}" |\
bcftools view \
	--threads ${threads} \
	--output-type z \
	--output ${out_vcf}
	
bcftools index --threads ${threads} --tbi ${out_vcf}
