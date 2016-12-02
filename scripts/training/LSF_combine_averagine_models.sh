#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_combine_averagine_models.sh
#BSUB -q day
#BSUB -n 1
#BSUB -M 1

module load python/2.7.1

source ../config.sh

python ${SOURCE_DIR}/scripts/training/combineModels.py ${SOURCE_DIR}/misc/schema.xml $AVERAGINE_OUT_DIR > ${AVERAGINE_OUT_DIR}"/models.xml"
