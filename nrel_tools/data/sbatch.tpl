#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --time={run_time}
#SBATCH --nodes=1
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err
#SBATCH --account={account}
{email}

# Load Gaussian module to set environment
module load gaussian/G16B
# Run gaussian Peregrine script (performs much of the Gaussian setup)
g16_eagle

cd $SLURM_SUBMIT_DIR

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss {job_name} {opt_old_name} -c {run_gauss_ini}

rm $SCRATCH2/*
rmdir $SCRATCH