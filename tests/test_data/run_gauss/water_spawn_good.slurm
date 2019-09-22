#!/bin/bash
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --job-name=water
#SBATCH --output=water.out
#SBATCH --error=water.err
#SBATCH --account=bpms
#SBATCH --qos=normal


# Load Gaussian module to set environment
module load gaussian/G16B
# Run gaussian Peregrine script (performs much of the Gaussian setup)
g16_eagle

cd $SLURM_SUBMIT_DIR

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss water -o water -c /Users/hmayes/bee/code/python/nrel_tools/water.ini

rm $SCRATCH2/*
rmdir $SCRATCH