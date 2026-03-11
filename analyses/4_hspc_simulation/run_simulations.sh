#!/bin/bash

#SBATCH -N 1
#SBATCH --time=6:00:00
#SBATCH --account=rmccoy22
#SBATCH --partition=parallel
#SBATCH --mem=1G
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-1000%50

# In each array job, run a simulation of HSPC population dynamics over the lifetime of one individual.

#========#
# Inputs #
#========#

# Path to SLiM binary
SLIM=~/code/SLiM-3.7.1/bin/slim
# Path to file of telomere lengths
TELOLENS=/scratch16/rmccoy22/syan11/clonal_hematopoesis/tertp_simulation/pop_telo_lens.txt
# SLiM simulation script
simfile=$1

# Directory and prefix for output files
outdir=$2
outpref=$3
outfile=$outdir/${outpref}_${mut_rate}_${SLURM_ARRAY_TASK_ID}.out

# Get nth line from file of telomere lengths to get starting TL for this simulation
line_num=${SLURM_ARRAY_TASK_ID}
telo_mean=$(sed -n "${line_num}p" $TELOLENS)

# Run simulation
$SLIM \
	-d mut_rate=$mut_rate \
	-d telo_mean=$telo_mean \
	-d shorten_mean=$shorten_mean \
	$simfile \
	> $outfile