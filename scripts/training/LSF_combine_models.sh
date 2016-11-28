#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_combine_models.sh
#BSUB -q day
#BSUB -n 1
#BSUB -M 1

module load python/2.7.1

source ../config.csh

python ${SOURCE_DIR}/scripts/training/combineModels.py ${SOURCE_DIR}/misc/schema.xml $ROOT_OUT_DIR > ${ROOT_OUT_DIR}"/models.xml"

