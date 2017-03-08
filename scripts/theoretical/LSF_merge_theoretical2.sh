#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_merge_theoretical2.sh
#BSUB -q week
#BSUB -o /netscr/dennisg/log/LSF_merge_theoretical.log.%J
#BSUB -n 1
#BSUB -M 1

module load python/2.7.1
module load r/3.2.2

source ../config.sh

set BIN_SIZE_RESIDUAL = 30
set BIN_SIZE_CHISQUARE = 30
set IN_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical/"

python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR "merged_precursor_scores_" F > ${IN_DIR}"/precursor_scores.txt"
python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR "merged_precursor_residuals_" F > ${IN_DIR}"/precursor_residuals.txt"

python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR "merged_fragment_scores_" T > ${IN_DIR}"/fragment_scores.txt"
python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR "merged_fragment_residuals_" T > ${IN_DIR}"/fragment_residuals.txt"

Rscript plotComparisons.R ${IN_DIR}"/precursor_residuals.txt" ${IN_DIR}"/precursor_residuals.eps" ${BIN_SIZE_RESIDUAL} F
Rscript plotComparisons.R ${IN_DIR}"/precursor_scores.txt" ${IN_DIR}"/precursor_chisquared.eps" ${BIN_SIZE_CHISQUARE} F

Rscript plotComparisons.R ${IN_DIR}"/fragment_residuals.txt" ${IN_DIR}"/fragment_residuals.eps" ${BIN_SIZE_RESIDUAL} T
Rscript plotComparisons.R ${IN_DIR}"/fragment_scores.txt" ${IN_DIR}"/fragment_chisquared.eps" ${BIN_SIZE_CHISQUARE} T
