#!/usr/bin/env bash
# Set script variables
INPUT_BASENAME=ethylrad
INPUT_FILE=tests/test_data/run_gauss/ethylrad.com
GAUSSIAN_EXEC=g16
MEMSIZE=5GB
SCRATCH=/scratch/hmayes/ethylrad_${SLURM_JOB_ID}
SCRATCH2=/dev/shm
INFILE=infile_${INPUT_BASENAME}

mkdir ${SCRATCH}
# Check on editing input file. If scratch directories
# are listed then file is used un-changed, if 3-line
# header not present, then script prepends these lines
# to the input file to be used in execution line
NUMRWFLINES=`grep "RWF" ${INPUT_FILE} | wc -l`
if [ ${NUMRWFLINES} -eq 1 ]; then
    echo "standard file found"
    cp ${INPUT_FILE} ${INFILE}
else
    echo "prepending lines to input file"
    echo "%RWF=${SCRATCH2}/,$MEMSIZE,${SCRATCH}/,-1" > ${INFILE}
    echo "%NoSave" >> ${INFILE}
    
    echo "%Chk=${SCRATCH2}/ethylrad.chk" >> ${INFILE}
    echo "%CPU=0-35" >> ${INFILE}
    echo "%Mem=72GB" >> ${INFILE}
    cat ${INPUT_FILE} >> ${INFILE}
fi


# Set required Gaussian environment variables
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then
    export GAUSS_LFLAGS='-vv -opt "Tsnet.Node.lindarsharg: ssh"'
    export GAUSS_EXEDIR=$g16root/g16/linda-exe:${GAUSS_EXEDIR}
fi
export GAUSS_SCRDIR=${SCRATCH2}

# Gaussian needs scratch directories
# If desired, make sure scratch is clear before starting
# rm ${SCRATCH2}/*

# Run Gaussian job
${GAUSSIAN_EXEC} < ${INFILE} >& ${INPUT_BASENAME}.log

rm ${INFILE}
cp ${SCRATCH2}/${INPUT_BASENAME}.chk .

# If desired, clean-up files or remove folder
# rm ${SCRATCH}/*
rm -r ${SCRATCH}