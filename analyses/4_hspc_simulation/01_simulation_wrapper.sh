#!/bin/bash

# Submits SLURM jobs to run the `run_simulations.sh` script, which
# runs a single simulation of HSPC population dynamics
# over the lifetime of one individual.

# `mut` is the per-division mutation rate at the TERTp site
for mut in 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9
do
	sbatch --job-name=popTL_${mut} \
		--export=ALL,mut_rate=$mut,shorten_mean=50 \
		run_simulations.sh \
		hsc_simulation_v10_shortTL_vaf.slim \
		output \
		popTL_unif
done