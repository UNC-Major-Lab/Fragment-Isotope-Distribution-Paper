#!/usr/bin/env bash

mkdir -p out/Average_Spline/spline/model/
mkdir out/Average_Spline/spline/gof/
mkdir out/Average_Spline/spline/hist/
mkdir out/Average_Spline/spline/res/
mkdir out/Average_Spline/spline/scatter/
mkdir out/Average_Spline/spline/eval/

for S in $(seq 0 5);
do
    set dir = "S$S"
    mkdir -p out/${dir}/spline/model/
    mkdir out/${dir}/spline/gof/
    mkdir out/${dir}/spline/hist/
    mkdir out/${dir}/spline/res/
    mkdir out/${dir}/spline/scatter/
    mkdir out/${dir}/spline/eval/
done

