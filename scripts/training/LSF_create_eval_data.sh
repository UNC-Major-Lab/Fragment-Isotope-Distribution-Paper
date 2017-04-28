#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_create_eval_data.sh
#BSUB -q hour
#BSUB -o /netscr/dennisg/log/create_eval_data.log.%J
#BSUB -n 1

module load gcc/4.8.1

source ../config.sh

set DATA_DIR = ${SPLINE_OUT_DIR}"/"
set OUT_DIR = ${ROOT_OUT_DIR}"/training/"

mkdir -p ${OUT_DIR}/proteome/averagine

${BUILD_DIR}/GenerateTrainingData $FASTA ${OUT_DIR}/proteome/ $MAX_SAMPLED_MASS_SULFUR $MAX_ISOTOPE_DEPTH $MONO
