#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_plots.sh
#BSUB -q hour
#BSUB -o /netscr/dennisg/log/training_plots.log.%J
#BSUB -n 1

module load r/3.2.2

source ../config.sh

set dir = "S0"
set DATA_DIR = ${SPLINE_OUT_DIR}"/"
set OUT_DIR = ${ROOT_OUT_DIR}"/training/"

Rscript plotModel.R ${DATA_DIR}"spline_eval/" ${DATA_DIR} 0 $MAX_SULFUR ${OUT_DIR}"Precursor0_model.eps"