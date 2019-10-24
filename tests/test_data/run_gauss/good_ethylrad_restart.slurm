#!/bin/bash
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=ethylrad_restart
#SBATCH --output=ethylrad_restart.out
#SBATCH --error=ethylrad_restart.err
#SBATCH --account=nwu150
#SBATCH --qos=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=heather.mayes@nrel.gov

# Load Gaussian module to set environment
module load gaussian/G16B

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss ethylrad_restart  -c ethylrad_restart.ini

rm $SCRATCH2/*
rmdir $SCRATCH