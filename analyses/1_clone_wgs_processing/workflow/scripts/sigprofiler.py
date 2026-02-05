#!/usr/bin/env python

# Run unsupervised SigProfiler assignment on VCF files

#Import libraries
import sys

from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall

#Check input format
if len(sys.argv) != 4:
        print("ERROR: Incorrect syntax")
        print("USAGE: ./Run_SigProfiler [input_dir] [output_dir] [genome_build]")
        print("EXITING...")

#Load variables
sample_dir = sys.argv[1]
output_dir = sys.argv[2]
genome_build = sys.argv[3]

#Install reference genome
genInstall.install(genome_build)

#Run SigProfiler
Analyze.cosmic_fit(samples = sample_dir, 
        output = output_dir,
        input_type = "vcf",
        context_type = "96",
        genome_build="GRCh38",
        verbose=True,
        make_plots=True,
        exome=False)