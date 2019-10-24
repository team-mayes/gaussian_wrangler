#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=ethylrad_opt
#SBATCH --output=ethylrad_opt.out
#SBATCH --error=ethylrad_opt.err
#SBATCH --account=bpms
#SBATCH --qos=high


# Load Gaussian module to set environment
module load gaussian/G16B

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss ethylrad_opt -o ethylrad_opt -c tests/test_data/run_gauss/run_spawn_all_new0.ini

rm $SCRATCH2/*
rmdir $SCRATCH