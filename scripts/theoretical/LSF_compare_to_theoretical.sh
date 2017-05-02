#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_compare_to_theoretical.sh[1-100]
#BSUB -q day
#BSUB -o /netscr/dennisg/log/LSF_compare_to_theoretical.log.%J
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1
module load r/3.2.2

source ../config.sh

set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical"

set BIN_SIZE_RESIDUAL = 0.0025
set BIN_SIZE_CHISQUARE = 0.025

mkdir -p $OUT_DIR

${BUILD_DIR}/CompareToTheoretical $FASTA $LSB_JOBINDEX 100 0 ${OUT_DIR}"/residuals_precursor_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/scores_precursor_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/stats_precursor_"${LSB_JOBINDEX}".out" $BIN_SIZE_CHISQUARE $BIN_SIZE_RESIDUAL
#${BUILD_DIR}/CompareToTheoretical $FASTA $LSB_JOBINDEX 100 1 ${OUT_DIR}"/residuals_fragment_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/scores_fragment_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/stats_fragment_"${LSB_JOBINDEX}".out" $BIN_SIZE_CHISQUARE $BIN_SIZE_RESIDUAL
