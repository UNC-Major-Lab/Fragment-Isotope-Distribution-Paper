#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J compare_to_theoretical.sh[1-300]
#BSUB -q hour
#BSUB -o /netscr/dennisg/log/compare_to_theoretical.log.%J
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1
module load r/3.2.2

source config.csh

set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical"

mkdir -p $OUT_DIR

${BUILD_DIR}/CompareToTheoretical $FASTA $LSB_JOBINDEX 300 > ${ROOT_OUT_DIR}"/compare_to_theoretical/"${LSB_JOBINDEX}".out"

