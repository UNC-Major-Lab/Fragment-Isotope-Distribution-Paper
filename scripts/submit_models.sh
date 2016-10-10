#!/bin/csh 

set S = $1
set CS = $2
set Se = $3
set CSe = $4 
set max_mass = $5

set ROOT_DIR = "/netscr/dennisg/MSFragmentDeisotoper/S${S}_CS${CS}_Se${Se}_CSe${CSe}"

foreach file (${ROOT_DIR}/data/*)
    if (-d $file) then
        echo "Skipping $file (is a directory)"
    else
	    bsub -q day -J "single.sh" -o "/netscr/dennisg/single.log.%J" -M 20 "./single_model.sh $file $S $CS $Se $CSe $max_mass"
    endif
end

