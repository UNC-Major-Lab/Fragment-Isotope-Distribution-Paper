#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J Isotopes.sh[1-16]
#BSUB -q day
#BSUB -o /netscr/dennisg/Isotopes.log.%J
#BSUB -n 1
#BSUB -M 20

module load gcc/4.8.1
module load r/3.2.2
module load matlab

source config.csh

echo $LSB_JOBINDEX

set index = $LSB_JOBINDEX
set CS = 0
set S = 0
set CSe = 0
set Se = 0

if ($LSB_JOBINDEX <= 1) then

else if ($LSB_JOBINDEX <= 2) then
    set S = 1
else if ($LSB_JOBINDEX <= 4) then
    set index = `expr $LSB_JOBINDEX - 2`
    set CS = `expr $index - 1`
    set S = `expr 3 - $index`
else if ($LSB_JOBINDEX <= 6) then
    set index = `expr $LSB_JOBINDEX - 4`
    set CS = `expr $index - 1`
    set S = `expr 4 - $index`
else if ($LSB_JOBINDEX <= 9) then
    set index = `expr $LSB_JOBINDEX - 6`
    set CS = `expr $index - 1`
    set S = `expr 5 - $index`
else if ($LSB_JOBINDEX <= 12) then
    set index = `expr $LSB_JOBINDEX - 9`
    set CS = `expr $index - 1`
    set S = `expr 6 - $index`
else if ($LSB_JOBINDEX <= 16) then
    set index = `expr $LSB_JOBINDEX - 12`
    set CS = `expr $index - 1`
    set S = `expr 7 - $index`
endif

set OUT_DIR = ${ROOT_OUT_DIR}"/S${S}_CS${CS}_Se${Se}_CSe${CSe}"

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

${BUILD_DIR}/PeptideFragmentSampler ${OUT_DIR}/data/ $MAX_SAMPLED_MASS $NUM_SAMPLES $S $CS $Se $CSe $MAX_ISOTOPE

chmod 775 ${OUT_DIR}/data/*

${SOURCE_DIR}/scripts/training/submit_models.sh $S $CS $Se $CSe $MAX_SAMPLED_MASS
