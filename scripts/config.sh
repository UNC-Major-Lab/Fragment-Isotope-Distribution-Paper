#!/bin/csh
set SOURCE_DIR="~/Fragment-Isotope-Distribution-Paper/"
set BUILD_DIR=${SOURCE_DIR}"/build/"
set DATA_DIR=${SOURCE_DIR}"/data/"

set ROOT_OUT_DIR="/netscr/dennisg/MSFragmentDeisotoper/"
set SULFUR_OUT_DIR=${ROOT_OUT_DIR}"/SulfurCorrected/"
set AVERAGINE_OUT_DIR=${ROOT_OUT_DIR}"/Averagine/"

set LOG_DIR="/netscr/dennisg/log/"

set FASTA=${DATA_DIR}"/human_sp_112816.fasta"

set NUM_SAMPLES="1000"
set MAX_SAMPLED_MASS="10000"
set MAX_ISOTOPE="20"
set MIN_PERCENTAGE_AVERAGINE="0.001"
