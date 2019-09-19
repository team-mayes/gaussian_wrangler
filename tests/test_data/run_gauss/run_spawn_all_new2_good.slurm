#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=ethylrad_opt
#SBATCH --output=ethylrad_opt.out
#SBATCH --error=ethylrad_opt.err
#SBATCH --account=bpms


# Load Gaussian module to set environment
module load gaussian/G16B
# Run gaussian Peregrine script (performs much of the Gaussian setup)
g16_eagle

cd $SLURM_SUBMIT_DIR

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss ethylrad_opt -o ethylrad_opt -c /Users/hmayes/bee/code/python/nrel_tools/tests/test_data/run_gauss/run_spawn_all_new2.ini

rm $SCRATCH2/*
rmdir $SCRATCH