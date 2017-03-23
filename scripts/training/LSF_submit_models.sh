#!/bin/csh 

source ../config.sh

set S = $1
if ($S < 0) then
    set dir = "Average_Spline"
else
    set $dir = "S$S"
endif

set max_mass = $2

set OUT_DIR = ${SPLINE_OUT_DIR}"/${dir}"

foreach file (${OUT_DIR}/data/*)
    if (-d $file) then
        echo "Skipping $file (is a directory)"
    else
	    bsub -q week -J "single.sh" -o ${LOG_DIR}"/single.log.%J" -M 20 ${SOURCE_DIR}"scripts/training/single_model.sh $file $S $max_mass $OUT_DIR"
    endif
end

