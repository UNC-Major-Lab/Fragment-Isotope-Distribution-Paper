#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_compare_to_targeted.sh
#BSUB -q day
#BSUB -o /netscr/dennisg/log/LSF_compare_to_targeted.log.%J
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1
module load r/3.2.2

source ../config.sh

set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/LT"

mkdir -p $OUT_DIR

${BUILD_DIR}/CompareToTargeted ${DATA_DIR}"/Neuro_04.mzML" ${DATA_DIR}"/Neuro_04_centroid.mzML" ${OUT_DIR}"/out04.tab" ${OUT_DIR}"/calc_out04.tab" ${OUT_DIR}"/scores_out04.tab"
Rscript lowThroughput/IndividualSpectrumIsotopes.R ${OUT_DIR}"/out04.tab" ${OUT_DIR}"/calc_out04.tab" ${OUT_DIR}"/scores_out04.tab" ${OUT_DIR}"/low_throughput.eps"