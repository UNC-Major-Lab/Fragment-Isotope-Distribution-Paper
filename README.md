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
```ShellSession
$ ./GenerateTrainingData 
 USAGE: GenerateTrainingData fasta_path out_path max_mass S num_samples max_depth mono
 
 fasta_path: The path to the fasta file to train the splines on.
 out_path: The path to the directory that will store the training data, e.g. ~/data/
 max_mass: maximum mass allowed for sampled peptides, e.g. 8500
 max_depth: The number of isotopes to generate training data for, e.g. 3 = M0,M1,M2
 mono: should monoisotopic masses be used or average? 1=mono, 0=average
 S: number of sulfurs that should be in the fragment ion. Use -1 for all (e.g. 0,1,2..)
 num_samples: number of random peptides to make for each peptide length
 
$ ./GenerateTrainingData data/human_sp_112816.fasta out/ 10000 5 1 Average_Spline 300
```

### Generate splines
Open MATLAB
Navigate to scripts/training folder

Call the IsotopeSpline function:

Usage: IsotopeSpline(knot_spacing, num_sulfurs, isotope, path_to_training_data, out_path_histogram, out_path_scatter_plot, out_path_residual_plot, out_path_GOF_stats, out_path_xml, out_path_spline_eval_for_figure1)  

To create the model for the monoisotope and peptides containing any number of sulfurs:
```Matlab
IsotopeSpline(1000,'Average_spline','0','out/Precursor0.txt','out/hist.pdf','out/scatter.pdf','out/res.pdf','out/gof.txt','out/model.xml','out/eval.tab')
```
To create the model for the M+1 isotope and peptides containing any number of sulfurs:
```Matlab
IsotopeSpline(1000,'Average_spline','1','out/Precursor0.txt','out/hist.pdf','out/scatter.pdf','out/res.pdf','out/gof.txt','out/model.xml','out/eval.tab')
```
To create the model for the monoisotope and peptides containing 0 sulfurs:
```Matlab
IsotopeSpline(1000,'0','0','out/Precursor0.txt','out/hist.pdf','out/scatter.pdf','out/res.pdf','out/gof.txt','out/model.xml','out/eval.tab')
```
The previous commands create separate spline model .xml files. Use this command to merge them into a singe .xml file.

Usage: combineModels.py path_to_spline_xmls max_isotope_depth max_sulfur
```ShellSession
$ python ../scripts/training/combineModels.py out/ 5 5 > out/IsotopeSplines.xml
```

### Figures S-1 and S-2

USAGE: plotModel.R path_to_eval_data data_dir isotope num_sulfur out_path_for_figure max_mass

Plot average spline for the monoisotope
```ShellSession
$ Rscript ../scripts/Training/plotModel.R ${DATA_DIR}"spline_eval/" ${DATA_DIR} 0 -1 Average_precursor0_model.eps 10000
```
Plot average spline for the M+1 isotope
```ShellSession
$ Rscript ../scripts/Training/plotModel.R ${DATA_DIR}"spline_eval/" ${DATA_DIR} 1 -1 Average_precursor1_model.eps 10000
```
Plot sulfur-specific spline for the monoisotope and 0 sulfurs
```ShellSession
$ Rscript ../scripts/Training/plotModel.R ${DATA_DIR}"spline_eval/" ${DATA_DIR} 0 0 precursor0_model.eps 10000
```

### Figure S-3

USAGE: SpeedTest min_mass max_mass max_sulfurs num_tests
```ShellSession
$ ./SpeedTest 400 9500 5 1e5 > out/runtimes.out
$ Rscript ../scripts/theoretical/plotRuntimeComparisons.R out/runtimes.out out/runtimes.eps
```

### Figure 1

### Figure 2
#${BUILD_DIR}/CompareToTheoretical $FASTA $LSB_JOBINDEX 100 1 ${OUT_DIR}"/residuals_fragment_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/scores_fragment_"${LSB_JOBINDEX}".out" ${OUT_DIR}"/stats_fragment_"${LSB_JOBINDEX}".out" $BIN_SIZE_CHISQUARE $BIN_SIZE_RESIDUAL
#Rscript plotComparisons.R ${IN_DIR}"/fragment_scores.txt" ${IN_DIR}"/fragment_chisquared.eps" ${BIN_SIZE_CHISQUARE} T


```ShellSession
$ ./CompareToTheoretical data/human_sp_112816.fasta 1 100 1 out/residuals_fragment".out" out/scores_fragment.out" out/stats_fragment.out" 0.1 0.0025
$ Rscript plotComparisons.R out/fragment_scores.txt out/fragment_chisquared.eps 0.1 T
```

### Figure 3

#${BUILD_DIR}/CompareToTargeted ${DATA_DIR}"/Neuro_04.mzML" ${DATA_DIR}"/Neuro_04_centroid.mzML" ${OUT_DIR}"/out04.tab" ${OUT_DIR}"/calc_out04.tab" ${OUT_DIR}"/scores_out04.tab"
#Rscript IndividualSpectrumIsotopes.R ${OUT_DIR}"/out04.tab" ${OUT_DIR}"/calc_out04.tab" ${OUT_DIR}"/scores_out04.tab" ${OUT_DIR}"/low_throughput.eps"


```ShellSession
$ ./CompareToTargeted data/Neuro_04.mzML data/Neuro_04_centroid.mzML out/out04.tab out/calc_out04.tab out/scores_out04.tab
$ Rscript IndividualSpectrumIsotopes.R out/out04.tab out/calc_out04.tab out/scores_out04.tab out/low_throughput.eps
```

### Figure 4 and Table 1

#CompareToShotgun ${DATA_DIR}"/HELA_2017-10-25_CID25_OT.mzML" ${DATA_DIR}"/HELA_2017-10-25_CID25_OT.idXML" 0.0 $OUT_DIR MS2 MS2 CID_25
```ShellSession
$ ./CompareToShotgun data/HELA_2017-10-25_CID25_OT.mzML data/HELA_2017-10-25_CID25_OT.idXML 0.0 out/ MS2 MS2 CID_25
```


