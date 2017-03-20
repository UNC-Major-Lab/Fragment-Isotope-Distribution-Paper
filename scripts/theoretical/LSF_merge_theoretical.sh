#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_merge_theoretical.sh[1-100]
#BSUB -q week
#BSUB -o /netscr/dennisg/log/LSF_merge_theoretical.log.%J
#BSUB -n 1
#BSUB -M 1

module load python/2.7.1 
module load r/3.2.2

source ../config.sh

set BIN_SIZE_RESIDUAL = 0.0025
set BIN_SIZE_CHISQUARE = 0.025
set IN_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical/"

python ${SOURCE_DIR}/scripts/theoretical/mergeTheoreticalResiduals.py $IN_DIR $BIN_SIZE_CHISQUARE "scores_precursor_" $LSB_JOBINDEX 100 > ${IN_DIR}"/merged_precursor_scores_"${LSB_JOBINDEX}".out"
python ${SOURCE_DIR}/scripts/theoretical/mergeTheoreticalResiduals.py $IN_DIR $BIN_SIZE_RESIDUAL "residuals_precursor_" $LSB_JOBINDEX 100 > ${IN_DIR}"/merged_precursor_residuals_"${LSB_JOBINDEX}".out"

python ${SOURCE_DIR}/scripts/theoretical/mergeTheoreticalComparisons.py $IN_DIR $BIN_SIZE_CHISQUARE "scores_fragment_" $LSB_JOBINDEX 100 > ${IN_DIR}"/merged_fragment_scores_"${LSB_JOBINDEX}".out"
python ${SOURCE_DIR}/scripts/theoretical/mergeTheoreticalComparisons.py $IN_DIR $BIN_SIZE_RESIDUAL "residuals_fragment_" $LSB_JOBINDEX 100 > ${IN_DIR}"/merged_fragment_residuals_"${LSB_JOBINDEX}".out"

