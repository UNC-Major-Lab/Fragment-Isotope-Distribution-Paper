#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_create_average_training_data.sh
#BSUB -q day
#BSUB -o /netscr/dennisg/Isotopes.log.%J
#BSUB -n 1
#BSUB -M 20

module load gcc/4.8.1
module load r/3.2.2
module load matlab

source ../config.csh

set OUT_DIR = ${ROOT_OUT_DIR}"/Average_Spline"

mkdir -p $OUT_DIR

rm -r ${OUT_DIR}/data/

mkdir ${OUT_DIR}/data/

rm -r ${OUT_DIR}/spline/

mkdir ${OUT_DIR}/spline/
mkdir ${OUT_DIR}/spline/3D/
mkdir ${OUT_DIR}/spline/model/
mkdir ${OUT_DIR}/spline/gof/
mkdir ${OUT_DIR}/spline/hist/
mkdir ${OUT_DIR}/spline/res3D/
mkdir ${OUT_DIR}/spline/scatter/

${BUILD_DIR}/PeptideFragmentSampler ${ROOT_OUT_DIR}/sulfur_distribution.tab ${OUT_DIR}/data/ $MAX_SAMPLED_MASS $MAX_ISOTOPE

chmod 775 ${OUT_DIR}/data/*
