#!/bin/bash
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --job-name=water
#SBATCH --output=water.out
#SBATCH --error=water.err
#SBATCH --account=nwu150
#SBATCH --qos=normal

#SBATCH --signal=B:USR1@30

copy_chk_before_exit()
{
    scp ${SLURM_JOB_NODELIST}:/dev/shm/*chk ${SLURM_SUBMIT_DIR}
    echo "function copy_chk_before_exit called at $(date)"
}
trap 'copy_chk_before_exit' EXIT
trap 'copy_chk_before_exit' USR1

# Load Gaussian module to set environment
module load gaussian/G16B

cd ${SLURM_SUBMIT_DIR}

run_gauss water -o water -c water.ini
