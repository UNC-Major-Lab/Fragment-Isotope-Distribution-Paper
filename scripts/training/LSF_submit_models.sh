#!/bin/csh 

source ../config.csh

set S = $1
set CS = $2
set Se = $3
set CSe = $4 
set max_mass = $5

set OUT_DIR = ${ROOT_OUT_DIR}"/S${S}_CS${CS}_Se${Se}_CSe${CSe}"

foreach file (${OUT_DIR}/data/*)
    if (-d $file) then
        echo "Skipping $file (is a directory)"
    else
	    bsub -q day -J "single.sh" -o ${LOG_DIR}"/single.log.%J" -M 20 ${SOURCE_DIR}"scripts/training/single_model.sh $file $S $CS $Se $CSe $max_mass $OUT_DIR"
    endif
end

