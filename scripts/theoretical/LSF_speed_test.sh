#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_speed_test.sh
#BSUB -q week
#BSUB -o /netscr/dennisg/log/LSF_speed_test.log.%J
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1
module load r/3.2.2

source ../config.sh

set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical/"

set MIN_MASS = 400
set MAX_MASS = 9500
set MAX_DEPTH = 11
set NUM_TESTS = 1e5

${BUILD_DIR}/SpeedTest $MIN_MASS $MAX_MASS $MAX_DEPTH $NUM_TESTS > ${OUT_DIR}"/runtimes.out"
Rscript plotRuntimeComparisons.R ${OUT_DIR}"/runtimes.out" ${OUT_DIR}"/runtimes.eps"