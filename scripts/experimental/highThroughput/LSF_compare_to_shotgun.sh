#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J LSF_compare_to_shotgun.sh
#BSUB -q day
#BSUB -o /netscr/dennisg/log/LSF_compare_to_shotgun.log.%J
#BSUB -n 1
#BSUB -M 8

module load gcc/4.8.1
module load r/3.2.2

source ../../config.sh

#################
# MS2 comparisons
########################################################################################################################
#set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_16_QuadIso_OT_Full_MS2"
#mkdir -p $OUT_DIR

#${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-01-31_5.mzML" ${DATA_DIR}"/HELA_2017-01-31_5.idXML" 0.0 $OUT_DIR "MS2"
#Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################

########################################################################################################################
#set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_16_QuadIso_IT_Full_MS2"
#mkdir -p $OUT_DIR

#${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-05-16_205_IT.mzML" ${DATA_DIR}"/HELA_2017-05-16_205_IT.idXML" 0.0 $OUT_DIR "MS2"
#Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################

########################################################################################################################
#set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_16_ITIso_OT_Full_MS2"
#mkdir -p $OUT_DIR

#${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-05-19_250_ITIso.mzML" ${DATA_DIR}"/HELA_2017-05-19_250_ITIso.idXML" 0.0 $OUT_DIR "MS2"
#Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################

########################################################################################################################
#set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_60_ITIso_OT_Full_MS2"
#mkdir -p $OUT_DIR

#${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-05-19_246_FullWidth_ITIso.mzML" ${DATA_DIR}"/HELA_2017-05-19_246_FullWidth_ITIso.idXML" 1.0 $OUT_DIR "MS2"
#Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################

########################################################################################################################
#set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_60_QuadIso_OT_Full_MS2"
#mkdir -p $OUT_DIR

#${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-05-18_242_FullWidth_QuadIso.mzML" ${DATA_DIR}"/HELA_2017-05-18_242_FullWidth_QuadIso.idXML" 1.0 $OUT_DIR "MS2"
#Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################



#########################
# Alternating Comparisons
########################################################################################################################
#set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_16_AlternatingIso_OT_Full_MS2"
#mkdir -p $OUT_DIR

#${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-05-22_296_alternating_iso.mzML" ${DATA_DIR}"/HELA_2017-05-22_296_alternating_iso.idXML" 0.0 $OUT_DIR "alternatingMS2Iso"
#Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################

########################################################################################################################
#set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_16_QuadIso_OT_AlternatingFullSim_MS2"
#mkdir -p $OUT_DIR

#${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-05-15_194_SIM_alternating.mzML" ${DATA_DIR}"/HELA_2017-05-15_194_SIM_alternating.idXML" 0.0 $OUT_DIR "alternatingMS2SIM" 400 900
#Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################

########################################################################################################################
set OUT_DIR = ${ROOT_OUT_DIR}"/compare_to_experimental/HT_HELA_16_QuadIso_OT_AlternatingFrag_MS2"
mkdir -p $OUT_DIR

${BUILD_DIR}/CompareToShotgun ${DATA_DIR}"/HELA_2017-06-05_62_alternating_frag.mzML" ${DATA_DIR}"/HELA_2017-06-05_62_alternating_frag.idXML" 0.0 $OUT_DIR "alternating MS2 HCD_vs_CID" 400 900
Rscript PlotShotgunResults.R ${OUT_DIR}"/distributionScores.out" $OUT_DIR
########################################################################################################################



