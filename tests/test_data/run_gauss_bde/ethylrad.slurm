#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=ethylrad
#SBATCH --output=std.out
#SBATCH --error=std.err
#SBATCH --account=bpms

# Load Gaussian module to set environment
module load gaussian/G16B
cd $SLURM_SUBMIT_DIR

# Set script variables
run_gauss_bde


rm $SCRATCH/*
rmdir $SCRATCH