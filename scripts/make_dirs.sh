#!/usr/bin/env bash

mkdir -p out/Average_Spline/spline/model/
mkdir out/Average_Spline/out/spline/gof/
mkdir out/Average_Spline/out/spline/hist/
mkdir out/Average_Spline/out/spline/res/
mkdir out/Average_Spline/out/spline/scatter/
mkdir out/Average_Spline/out/spline/eval/

for S in $(seq 0 5);
do
    set dir = "S$S"
    mkdir -p out/${dir}/spline/model/
    mkdir out/${dir}/out/spline/gof/
    mkdir out/${dir}/out/spline/hist/
    mkdir out/${dir}/out/spline/res/
    mkdir out/${dir}/out/spline/scatter/
    mkdir out/${dir}/out/spline/eval/
done

