#!/usr/bin/env bash
#SBATCH --job-name=cellphy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=72:00:00

source ~/.bashrc
source activate CellPhy

~/bin/cellphy/cellphy.sh \
    FULL \
    -a \
    -p cellphy \
    /path/to/input.vcf
