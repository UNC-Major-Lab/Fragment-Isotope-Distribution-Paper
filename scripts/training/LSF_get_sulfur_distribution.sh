#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_get_sulfur_distribution.sh
#BSUB -q week
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1

source ../config.sh

mkdir -p ${SPLINE_OUT_DIR}

${BUILD_DIR}/GetSulfurDistribution $FASTA > ${SPLINE_OUT_DIR}"/sulfur_distribution.tab"

