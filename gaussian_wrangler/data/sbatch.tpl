#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --time={run_time}
#SBATCH --nodes=1
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err
#SBATCH --account={account}
#SBATCH --qos={qos}
{email}
#SBATCH --signal=B:USR1@30

copy_chk_before_exit()
{
    scp ${SLURM_JOB_NODELIST}:/dev/shm/*chk $SLURM_SUBMIT_DIR
    echo "function copy_chk_before_exit called at $(date)"
}
trap 'copy_chk_before_exit' EXIT
trap 'copy_chk_before_exit' USR1

# Load Gaussian module to set environment
module load gaussian/G16B

SCRATCH=/tmp/scratch/$SLURM_JOB_ID
mkdir $SCRATCH

cd $SLURM_SUBMIT_DIR

run_gauss {job_name} {opt_old_name} -c {run_gauss_ini}

rm $SCRATCH2/*
rmdir $SCRATCH