#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --job-name=ethylrad_opt_freq
#SBATCH --output=ethylrad_opt_freq.out
#SBATCH --error=ethylrad_opt_freq.err
#SBATCH --account=bpms
#SBATCH --qos=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=heather.mayes@nrel.gov
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

run_gauss ethylrad_opt  -c ethylrad_opt_freq.ini
