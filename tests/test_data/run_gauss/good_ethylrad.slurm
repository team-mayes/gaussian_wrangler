#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=ethylrad
#SBATCH --output=ethylrad.out
#SBATCH --error=ethylrad.err
#SBATCH --account=bpms
#SBATCH --qos=high


# Load Gaussian module to set environment
module load gaussian/G16B

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss ethylrad  -c ethylrad.ini

rm $SCRATCH2/*
rmdir $SCRATCH