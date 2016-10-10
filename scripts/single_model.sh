#!/bin/csh

set file = $1
set S = $2
set CS = $3
set Se = $4
set CSe = $5
set max_mass = $6

set ROOT_DIR = "/netscr/dennisg/MSFragmentDeisotoper/S${S}_CS${CS}_Se${Se}_CSe${CSe}"

set precursor = `echo "$file" | sed -n "s/.*\/[a-zA-Z]*\([0-9]*\)_[a-zA-Z]*\([0-9]*\).*/\1/p"`
set fragment = `echo "$file" | sed -n "s/.*\/[a-zA-Z]*\([0-9]*\)_[a-zA-Z]*\([0-9]*\).*/\2/p"`

set out_spline = ${ROOT_DIR}/spline/3D/Precursor${precursor}_Fragment${fragment}.pdf
set out_res_hist = ${ROOT_DIR}/spline/hist/Precursor${precursor}_Fragment${fragment}.pdf
set out_scatter = ${ROOT_DIR}/spline/scatter/Precursor${precursor}_Fragment${fragment}.tiff
set out_res3D = ${ROOT_DIR}/spline/res3D/Precursor${precursor}_Fragment${fragment}.tiff
set out_gof = ${ROOT_DIR}/spline/gof/Precursor${precursor}_Fragment${fragment}.txt
set out_models = ${ROOT_DIR}/spline/model/Precursor${precursor}_Fragment${fragment}.xml

set param = "FragmentIsotopeTensorProductSpline('"${max_mass}"','"${S}"','"${CS}"','"${Se}"','"${CSe}"','"${precursor}"','"${fragment}"','"${file}"','"${out_spline}"','"${out_res_hist}"','"${out_scatter}"','"${out_res3D}"','"${out_gof}"','"${out_models}"');quit force"

matlab -nodisplay -nosplash -nodesktop -r $param
