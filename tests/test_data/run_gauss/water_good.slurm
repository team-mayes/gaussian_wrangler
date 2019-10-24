#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=water
#SBATCH --output=water.out
#SBATCH --error=water.err
#SBATCH --account=bpms
#SBATCH --qos=high


# Load Gaussian module to set environment
module load gaussian/G16B

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss water  -c water.ini

rm $SCRATCH2/*
rmdir $SCRATCH