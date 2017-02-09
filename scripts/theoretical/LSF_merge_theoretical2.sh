#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_merge_theoretical.sh
#BSUB -q week
#BSUB -o /netscr/dennisg/log/LSF_merge_theoretical.log.%J
#BSUB -n 1
#BSUB -M 1

module load python/2.7.1
module load r/3.2.2

source ../config.sh

set BIN_SIZE = 0.0025
set IN_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical/"

python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR $BIN_SIZE "merged_precursor_scores_" 1 10 > ${IN_DIR}"/precursor_scores.txt"
python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR $BIN_SIZE "merged_precursor_residuals_" 1 10 > ${IN_DIR}"/precursor_residuals.txt"

python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR $BIN_SIZE "merged_fragment_scores_" 1 10 > ${IN_DIR}"/fragment_scores.txt"
python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR $BIN_SIZE "merged_fragment_residuals_" 1 10 > ${IN_DIR}"/fragment_residuals.txt"

Rscript plotComparisons.R ${IN_DIR}"/precursor_residuals.txt" ${IN_DIR}"/residuals.pdf"

