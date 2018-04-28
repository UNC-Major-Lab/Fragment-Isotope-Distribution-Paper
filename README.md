# Fragment-Isotope-Distribution-Paper
Source code associated with the paper 'Approximating isotopic distributions of biomolecule fragments'

This repository contains the source code necessary to reproduce the results from the paper.
It consists of the code to generate the training data for splines, generate the spline models,
and to compare approximate isotopic distributions to various theoretical and experimental isotopic distributions.

The scripts in the scripts/ directory are tailored to our LSF cluster to perform parallel job execution.
It's not neccessary to use these scripts, but they contain examples of how to execute the programs
and the values we used for various parameters.

## Requirements
1. CMake 2.8.3+
2. MATLAB R2016a
3. R 3.2.2
4. Python 2.7.1
5. OpenMS 2.1 (our custom fork, that includes the spline models. see below for instructions)

## Build Instructions
1. Build our version of OpenMS 2.1

Follow the instructions to be build OpenMS here: https://github.com/OpenMS/OpenMS/wiki/Building-OpenMS
but clone our version of OpenMS instead of theirs:

```ShellSession
$ git clone https://github.com/DennisGoldfarb/OpenMS.git
$ cd OpenMS
$ git checkout IsotopesSplines
```

2. Build this project

```ShellSession
$ git clone https://github.com/UNC-Major-Lab/Fragment-Isotope-Distribution-Paper.git
$ cd Fragment-Isotope-Distribution-Paper
$ mkdir build
$ cd build
$ cmake ../ -DOpenMS_DIR=/PATH/TO/OPENMS/BUILD
$ make
```

## Execution

### Generate Training data

#GenerateTrainingData FASTA_FILE OUT_PATH MAX_SAMPLED_MASS MAX_ISOTOPE_DEPTH MONO_MASS? NUM_SULFURS NUM_SAMPLES
```ShellSession
GenerateTrainingData data/human_sp_112816.fasta out/ 10000 5 1 Average_Spline 300
```

### Generate splines
Open MATLAB
Navigate to scripts/training folder
#IsotopeSpline("${SPLINE_BREAKS_SIZE}",'"${S}"','"${precursor}"','"${file}"','"${out_res_hist}"','"${out_scatter}"','"${out_res}"','"${out_gof}"','"${out_models}"','"${out_spline_eval}"')
IsotopeSpline(1000,'Average_spline','0','out/Precursor0.txt','out/hist.pdf','out/scatter.pdf','out/res.pdf','out/gof.txt','out/model.xml','out/eval.tab')

#python combineModels.py $SPLINE_OUT_DIR $MAX_ISOTOPE_DEPTH $MAX_SULFUR > ${SPLINE_OUT_DIR}"/IsotopeSplines.xml"
```ShellSession
python combineModels.py out/ 5 5 > out/IsotopeSplines.xml
```

### Figure 1
#GenerateTrainingData $FASTA ${OUT_DIR}/proteome/ $MAX_SAMPLED_MASS_SULFUR $MAX_ISOTOPE_DEPTH $MONO
#Rscript plotModel.R ${DATA_DIR}"spline_eval/" ${DATA_DIR} $PRECURSOR_ISOTOPE -1 ${OUT_DIR}"Average_precursor"${PRECURSOR_ISOTOPE}"_model.eps" 100000


```ShellSession
GenerateTrainingData data/human_sp_112816.fasta out/proteome/ 10000 5 1
Rscript plotModel.R ${DATA_DIR}"spline_eval/" ${DATA_DIR} $PRECURSOR_ISOTOPE -1 ${OUT_DIR}"Average_precursor"${PRECURSOR_ISOTOPE}"_model.eps" 100000
```

### Figure 2
#${BUILD_DIR}/CompareToTheoretical $FASTA $LSB_JOBINDEX 100 1 ${OUT_DIR}"/residuals_fragment_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/scores_fragment_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/stats_fragment_"${LSB_JOBINDEX}".out" $BIN_SIZE_CHISQUARE $BIN_SIZE_RESIDUAL
#Rscript plotComparisons.R ${IN_DIR}"/fragment_scores.txt" ${IN_DIR}"/fragment_chisquared.eps" ${BIN_SIZE_CHISQUARE} T


```ShellSession
CompareToTheoretical data/human_sp_112816.fasta 1 100 1 out/residuals_fragment".out" out/scores_fragment.out" out/stats_fragment.out" 0.1 0.0025
Rscript plotComparisons.R out/fragment_scores.txt out/fragment_chisquared.eps 0.1 T
```

### Figure 3

#${BUILD_DIR}/CompareToTargeted ${DATA_DIR}"/Neuro_04.mzML" ${DATA_DIR}"/Neuro_04_centroid.mzML" ${OUT_DIR}"/out04.tab" ${OUT_DIR}"/calc_out04.tab" ${OUT_DIR}"/scores_out04.tab"
#Rscript IndividualSpectrumIsotopes.R ${OUT_DIR}"/out04.tab" ${OUT_DIR}"/calc_out04.tab" ${OUT_DIR}"/scores_out04.tab" ${OUT_DIR}"/low_throughput.eps"


```ShellSession
${BUILD_DIR}/CompareToTargeted data/Neuro_04.mzML data/Neuro_04_centroid.mzML out/out04.tab out/calc_out04.tab out/scores_out04.tab
Rscript IndividualSpectrumIsotopes.R out/out04.tab out/calc_out04.tab out/scores_out04.tab out/low_throughput.eps
```

### Figure 4 and Table 1

#CompareToShotgun ${DATA_DIR}"/HELA_2017-10-25_CID25_OT.mzML" ${DATA_DIR}"/HELA_2017-10-25_CID25_OT.idXML" 0.0 $OUT_DIR MS2 MS2 CID_25
```ShellSession
CompareToShotgun data/HELA_2017-10-25_CID25_OT.mzML data/HELA_2017-10-25_CID25_OT.idXML 0.0 out/ MS2 MS2 CID_25
```

### Supplemental speed comparison

```ShellSession
SpeedTest 400 9500 5 1e5 > out/runtimes.out
Rscript plotRuntimeComparisons.R out/runtimes.out out/runtimes.eps
```