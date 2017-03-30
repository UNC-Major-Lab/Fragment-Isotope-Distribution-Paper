#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_create_training_data.sh[1-7]
#BSUB -q week
#BSUB -o /netscr/dennisg/log/Isotopes.log.%J
#BSUB -n 1

module load gcc/4.8.1
module load r/3.2.2
module load matlab

source ../config.sh

set S = `expr $LSB_JOBINDEX - 2`

if ($S < 0) then
    set dir = "Average_Spline"
else
    set dir = "S$S"
endif

set OUT_DIR = ${SPLINE_OUT_DIR}"/${dir}"

mkdir -p $OUT_DIR

rm -r ${OUT_DIR}/data/

mkdir ${OUT_DIR}/data/

rm -r ${OUT_DIR}/spline/

mkdir ${OUT_DIR}/spline/
mkdir ${OUT_DIR}/spline/model/
mkdir ${OUT_DIR}/spline/gof/
mkdir ${OUT_DIR}/spline/hist/
mkdir ${OUT_DIR}/spline/res/
mkdir ${OUT_DIR}/spline/scatter/

${BUILD_DIR}/GenerateTrainingData $FASTA ${OUT_DIR}/data/ $MAX_SAMPLED_MASS $S $MAX_ISOTOPE_DEPTH $MONO

chmod 775 ${OUT_DIR}/data/*

${SOURCE_DIR}/scripts/training/LSF_submit_models.sh $S $MAX_SAMPLED_MASS
