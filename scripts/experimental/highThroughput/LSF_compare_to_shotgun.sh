#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_compare_to_shotgun.sh
#BSUB -q week
#BSUB -o /netscr/dennisg/log/LSF_compare_to_shotgun.log.%J
#BSUB -n 1
#BSUB -M 8

module load gcc/4.8.1
module load r/3.2.2

source ../config.sh

set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT"

mkdir -p $OUT_DIR

${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-01-31_5.mzML" ${DATA_DIR}"/HELA_2017-01-31_5.idXML" 0.0 $OUT_DIR
Rscript highThroughput/PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR