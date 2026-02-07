#!/usr/bin/env bash

# This script will align paired-end reads to a reference genome using BWA mem, 
# compress the output, and index the compressed output. Readgroup information is either provided
# directly as command-line arguments or extracted from the fastq headers.

#--- LOAD RUN VARIABLES ---#

# command-line argument version
R1=$1 #Argument 1 = path to R1 fastq
R2=$2 #Argument 2 = path to R2 fastq
REFERENCE_PATH=$3 #Argument 3 = path to reference genome
OUTPUT=$4 #Argument 4 = path to output file

THREADS=$5 #Argument 5 = number of threads

# read group information
RGLB=$6 #Argument 6 = Library name
RGSM=$7 #Argument 7 = Sample name
RGPL=ILLUMINA #Platform 
RGID=`zcat ${R1} | head -n 1 | awk '{print $1}' | awk -F":" '{OFS="."; print $3,$4}' | xargs` #Read flowcell barcode and lane from first read header in $R1 to create ID

flowcell_barcode=`zcat ${R1} | head -n 1 | awk '{print $1}' | awk -F":" '{print $3}' | xargs` #Read flowcell barcode from $R1 for RGPU
lane=`zcat ${R1} | head -n 1 | awk '{print $1}' | awk -F":" '{print $4}' | xargs` #Read lane from $R1 for RGPU
sample_barcode=`zcat ${R1} | head -n 1 | awk '{print $2}' | awk -F":" '{print $(NF)}' | xargs` #Read sample barcode from $R1 for RGPU

#--- QC PRINTING ---#
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "REFERENCE PATH: ${REFERENCE_PATH}"
echo "THREADS: ${THREADS}"
echo "READ GROUP ID: ${RGID}"
echo "READ GROUP PLATFORM: ${RGPL}"
echo "READ GROUP PLATFORM UNIT: ${flowcell_barcode}.${lane}.${sample_barcode}"
echo "READ GROUP SAMPLE: ${RGSM}"
echo "READ GROUP LIBRARY: ${RGLB}"

#--- BWA ALIGNMENT, COMPRESSION, and INDEXING---#

# because we have nested variables, we need to use eval to expand the variables
CMD=`echo "bwa mem -t $THREADS -R $(echo "'@RG\tID:${RGID}\tPL:${RGPL}\tPU:${flowcell_barcode}.${lane}.${sample_barcode}\tLB:${RGLB}\tSM:${RGSM}'") ${REFERENCE_PATH} ${R1} ${R2}"`
echo "Running: ${CMD}"

eval ${CMD} | samtools view -Shb -@ $THREADS -o ${OUTPUT}

