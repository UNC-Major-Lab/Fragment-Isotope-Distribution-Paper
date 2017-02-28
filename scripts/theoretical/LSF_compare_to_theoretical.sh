#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_compare_to_theoretical.sh[1-100]
#BSUB -q week
#BSUB -o /netscr/dennisg/log/LSF_compare_to_theoretical.log.%J
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1
module load r/3.2.2

source ../config.sh

set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical"

mkdir -p $OUT_DIR

${BUILD_DIR}/CompareToTheoretical $FASTA $LSB_JOBINDEX 100 0 ${OUT_DIR}"/residuals_precursor_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/scores_precursor_"${LSB_JOBINDEX}".out"
${BUILD_DIR}/CompareToTheoretical $FASTA $LSB_JOBINDEX 100 1 ${OUT_DIR}"/residuals_fragment_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/scores_fragment_"${LSB_JOBINDEX}".out"

