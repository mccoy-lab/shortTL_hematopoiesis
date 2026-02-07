#!/usr/bin/env bash

# Run VQSR on VCF

while getopts ":m:a:d:t:h:o:x:i:c:p:v:" opt; do
	case ${opt} in
		m ) MILLS=$OPTARG
			;;
		a ) AXIOM=$OPTARG
			;;
		d ) DBSNP=$OPTARG
			;;
		t ) TGP=$OPTARG
			;;
		h ) HAPMAP=$OPTARG
			;;
		o ) OMNI=$OPTARG
			;;
		x ) MEMORY="-Xmx${OPTARG}m"
			;;
		i ) INPUT_VCF=$OPTARG
			;;
		c ) CHROMOSOME=$OPTARG
			;;
		p ) TMP_DIR=$OPTARG
			;;
		v ) OUTPUT_VCF=$OPTARG
			;;
		\? ) echo "Usage: cmd [-m] [-a] [-d] [-t] [-h] [-o] [-x memory] [-i input_vcf] [-c chromosome] [-p tmp_dir] [-v output_vcf]"
			exit 1
			;;
	esac
done

#Remove sample information
gatk MakeSitesOnlyVcf \
	-I ${INPUT_VCF} \
	-O ${TMP_DIR}/chr${CHROMOSOME}.sites_only.vcf.gz

#Indel filter
gatk --java-options ${MEMORY} VariantRecalibrator \
	-V ${TMP_DIR}/chr${CHROMOSOME}.sites_only.vcf.gz \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -an BaseQRankSum \
	-mode INDEL \
	--max-gaussians 4 \
	--resource:mills,known=false,training=true,truth=true,prior=12 ${MILLS} \
 	--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${AXIOM} \
 	--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${DBSNP} \
 	-O ${TMP_DIR}/chr${CHROMOSOME}.sites_only.indels.recal \
 	--tranches-file ${TMP_DIR}/chr${CHROMOSOME}.sites_only.indels.tranches

 #SNP filter
 gatk --java-options ${MEMORY} VariantRecalibrator \
	-V ${TMP_DIR}/chr${CHROMOSOME}.sites_only.vcf.gz \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
	-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -an BaseQRankSum \
	-mode SNP \
	--max-gaussians 6 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15 ${HAPMAP} \
    --resource:omni,known=false,training=true,truth=true,prior=12 ${OMNI} \
    --resource:1000G,known=false,training=true,truth=false,prior=10 ${TGP} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${DBSNP} \
 	-O ${TMP_DIR}/chr${CHROMOSOME}.sites_only.snps.recal \
 	--tranches-file ${TMP_DIR}/chr${CHROMOSOME}.sites_only.snps.tranches

if [[ ! -f ${TMP_DIR}/chr${CHROMOSOME}.sites_only.indels.recal ]]; then
	echo "ERROR: VariantRecalibrator Indel output not found."
	exit
fi

if [[ ! -f ${TMP_DIR}/chr${CHROMOSOME}.sites_only.snps.recal ]]; then
	echo "ERROR: VariantRecalibrator SNPs output not found."
	exit
fi

#Apply Indel VQSR
gatk --java-options ${MEMORY} ApplyVQSR \
	-V ${INPUT_VCF} \
	--recal-file ${TMP_DIR}/chr${CHROMOSOME}.sites_only.indels.recal \
	--tranches-file ${TMP_DIR}/chr${CHROMOSOME}.sites_only.indels.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	-mode INDEL \
	-O ${TMP_DIR}/chr${CHROMOSOME}.indel.recalibrated.vcf.gz

#Apply SNP VQSR
gatk --java-options ${MEMORY} ApplyVQSR \
	-V ${TMP_DIR}/chr${CHROMOSOME}.indel.recalibrated.vcf.gz \
	--recal-file ${TMP_DIR}/chr${CHROMOSOME}.sites_only.snps.recal \
	--tranches-file ${TMP_DIR}/chr${CHROMOSOME}.sites_only.snps.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	-mode SNP \
	-O ${OUTPUT_VCF}

# Remove intermediate files
rm ${TMP_DIR}/chr${CHROMOSOME}.sites_only.vcf.gz*
rm ${TMP_DIR}/chr${CHROMOSOME}.sites_only.indels.recal*
rm ${TMP_DIR}/chr${CHROMOSOME}.sites_only.indels.tranches
rm ${TMP_DIR}/chr${CHROMOSOME}.sites_only.snps.recal*
rm ${TMP_DIR}/chr${CHROMOSOME}.sites_only.snps.tranches
rm ${TMP_DIR}/chr${CHROMOSOME}.indel.recalibrated.vcf.gz*