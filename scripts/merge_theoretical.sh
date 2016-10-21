#!/bin/csh
#BSUB -L /bin/csh
#BSUB -J merge_theoretical.sh
#BSUB -q day
#BSUB -o /netscr/dennisg/log/merge_theoretical.log.%J
#BSUB -n 1
#BSUB -M 1

module load python/2.7.1 
module load r/3.2.2

source config.csh

set BIN_SIZE = 0.01
set IN_DIR = ${ROOT_OUT_DIR}"/compare_to_theoretical/"

python ${SOURCE_DIR}/scripts/mergeTheoreticalComparisons.py $IN_DIR $BIN_SIZE > ${IN_DIR}"/hist.txt"

