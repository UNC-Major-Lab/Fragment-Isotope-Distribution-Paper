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

set BIN_SIZE = 0.01
set IN_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical/"

python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR $BIN_SIZE "merged_scores_" 1 10 > ${IN_DIR}"/scores.txt"
python ${SOURCE_DIR}/scripts/theoretical/mergeHistogram.py $IN_DIR $BIN_SIZE "merged_residuals_" 1 10 > ${IN_DIR}"/residuals.txt"

Rscript plotComparisons.R ${IN_DIR}"/residuals.txt" ${IN_DIR}"/residuals.pdf"

