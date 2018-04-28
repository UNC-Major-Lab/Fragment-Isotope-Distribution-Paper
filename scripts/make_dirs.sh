#!/bin/csh

mkdir -p out/Average_Spline/spline/model/
mkdir -p out/Average_Spline/spline/gof/
mkdir -p out/Average_Spline/spline/hist/
mkdir -p out/Average_Spline/spline/res/
mkdir -p out/Average_Spline/spline/scatter/
mkdir -p out/Average_Spline/spline/eval/
mkdir -p out/Average_Spline/data/

mkdir -p out/proteome/
mkdir -p out/proteome/averagine/

set S = 0
while  ($S <= 5);
    set out_dir="S${S}"
    mkdir -p out/${out_dir}/spline/model/
    mkdir -p out/${out_dir}/spline/gof/
    mkdir -p out/${out_dir}/spline/hist/
    mkdir -p out/${out_dir}/spline/res/
    mkdir -p out/${out_dir}/spline/scatter/
    mkdir -p out/${out_dir}/spline/eval/
    mkdir -p out/${out_dir}/data/
    @ S++
end

