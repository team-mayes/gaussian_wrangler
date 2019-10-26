#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=water_opt_stable
#SBATCH --output=water_opt_stable.out
#SBATCH --error=water_opt_stable.err
#SBATCH --account=bpms
#SBATCH --qos=high

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

run_gauss water  -c water_opt_stable.ini
