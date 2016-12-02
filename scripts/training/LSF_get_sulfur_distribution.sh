#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_get_sulfur_distribution.sh
#BSUB -q hour
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1

source ../config.csh

${BUILD_DIR}/GetSulfurDistribution $FASTA > ${AVERAGINE_OUT_DIR}"/sulfur_distribution.tab"

