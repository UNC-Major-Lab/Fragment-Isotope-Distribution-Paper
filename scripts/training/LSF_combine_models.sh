#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_combine_models.sh
#BSUB -q week
#BSUB -o /netscr/dennisg/log/LSF_combine_models.log.%J
#BSUB -n 1
#BSUB -M 1

module load python/2.7.1

source ../config.sh

python ${SOURCE_DIR}/scripts/training/combineModels.py $SPLINE_OUT_DIR $MAX_ISOTOPE $MAX_SULFUR > ${SPLINE_OUT_DIR}"/IsotopeSplines.xml"