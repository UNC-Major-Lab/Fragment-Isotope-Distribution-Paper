#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_get_sulfur_distribution.sh
#BSUB -o /netscr/dennisg/log/LSF_get_sulfur_distribution.log.%J
#BSUB -q hour
#BSUB -n 1
#BSUB -M 1

module load gcc/4.8.1

source ../config.sh

mkdir -p ${SPLINE_OUT_DIR}

${BUILD_DIR}/GetSulfurDistribution $FASTA > ${SPLINE_OUT_DIR}"/sulfur_distribution.tab"

