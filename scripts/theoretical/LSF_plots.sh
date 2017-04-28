#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_plots.sh[1-21]
#BSUB -q hour
#BSUB -o /netscr/dennisg/log/theoretical_plots.log.%J
#BSUB -n 1

module load r/3.2.2

source ../config.sh

set DATA_DIR = ${SPLINE_OUT_DIR}"/"
set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical/"

set PRECURSOR_ISOTOPE = `expr $LSB_JOBINDEX - 1`

Rscript plotModelToProteome.R ${DATA_DIR}"spline_eval/" ${ROOT_OUT_DIR} $PRECURSOR_ISOTOPE $MAX_SULFUR ${OUT_DIR}"spline_comparison"${PRECURSOR_ISOTOPE}"_model.eps" 10000
