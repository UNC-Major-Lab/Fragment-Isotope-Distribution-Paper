#!/bin/csh

source ../config.sh

set max_mass = $1

set OUT_DIR = ${AVERAGINE_OUT_DIR}"/Average_Spline"

foreach file (${OUT_DIR}/data/*)
    if (-d $file) then
        echo "Skipping $file (is a directory)"
    else
	    bsub -q week -J "single.sh" -o ${LOG_DIR}"/single.log.%J" -M 20 ${SOURCE_DIR}"scripts/training/single_model.sh $file NA $max_mass $OUT_DIR"
    endif
end

