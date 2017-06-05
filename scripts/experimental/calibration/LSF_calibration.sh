#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_calibration.sh
#BSUB -q day
#BSUB -o /netscr/dennisg/log/LSF_calibration.log.%J
#BSUB -n 1
#BSUB -M 8

module load gcc/4.8.1
module load r/3.2.2

source ../../config.sh


########################################################################################################################
set OUT_DIR = ${ROOT_OUT_DIR}"/calibration/"
mkdir -p $OUT_DIR
########################################################################################################################


${BUILD_DIR}/ProcessCalibration ${DATA_DIR}"/OT_QuadIsolation.mzML" ${DATA_DIR}"/cal.txt" ${OUT_DIR}"OT_Quad.out" 200
Rscript PlotCalibration.R ${OUT_DIR}"/OT_Quad.out" ${OUT_DIR}"/OT_Quad.pdf"


${BUILD_DIR}/ProcessCalibration ${DATA_DIR}"/IT_QuadIsolation.mzML" ${DATA_DIR}"/cal.txt" ${OUT_DIR}"IT_Quad.out" 2000
Rscript PlotCalibration.R ${OUT_DIR}"/IT_Quad.out" ${OUT_DIR}"/IT_Quad.pdf"


${BUILD_DIR}/ProcessCalibration ${DATA_DIR}"/OT_ITIsolation.mzML" ${DATA_DIR}"/cal.txt" ${OUT_DIR}"OT_IT.out" 200
Rscript PlotCalibration.R ${OUT_DIR}"/OT_IT.out" ${OUT_DIR}"/OT_IT.pdf"


${BUILD_DIR}/ProcessCalibration ${DATA_DIR}"/IT_ITIsolation.mzML" ${DATA_DIR}"/cal.txt" ${OUT_DIR}"IT_IT.out" 2000
Rscript PlotCalibration.R ${OUT_DIR}"/IT_IT.out" ${OUT_DIR}"/IT_IT.pdf"


