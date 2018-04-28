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
6. gcc-4.8.1 and g++-4.8.1 (other version might work too, but that's what I used)

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
$ ../scripts/make_dirs.sh
$ ./GenerateTrainingData 
 USAGE: GenerateTrainingData fasta_path out_path max_mass max_depth mono S num_samples 
 or 
 GenerateTrainingData fasta_path out_path max_mass max_depth mono
 
 fasta_path: The path to the fasta file to train the splines on.
 out_path: The path to the directory that will store the training data, e.g. ~/data/
 max_mass: maximum mass allowed for sampled peptides, e.g. 8500
 max_depth: The number of isotopes to generate training data for, e.g. 3 = M0,M1,M2
 mono: should monoisotopic masses be used or average? 1=mono, 0=average
 S: number of sulfurs that should be in the fragment ion. Use -1 for all (e.g. 0,1,2..)
 num_samples: number of random peptides to make for each peptide length
 
$ ./GenerateTrainingData ../data/human_sp_112816.fasta out/Average_Spline/data/ 10000 5 1 -1 300
$ ./GenerateTrainingData ../data/human_sp_112816.fasta out/S0/data/ 10000 5 1 0 300
```

The above generates the training data for the first 5 isotopes for the average spline up to 10kDa and the sulfur-specific spline with 0 sulfurs. Repeat for the last command with different values of S for other sulfur-specific models. For publication we did the first 100 isotopes.

### Generate splines
Open MATLAB
Navigate to scripts/training folder

Call the IsotopeSpline function:

Usage: IsotopeSpline(knot_spacing, num_sulfurs, isotope, path_to_training_data, out_path_histogram, out_path_scatter_plot, out_path_residual_plot, out_path_GOF_stats, out_path_xml, out_path_spline_eval_for_figure1)  

To create the model for the monoisotope and peptides containing any number of sulfurs:
```Matlab
IsotopeSpline(1000,'Average_spline','0','../../build/out/Average_Spline/data/Precursor0.tab','../../build/out/Average_Spline/spline/hist/Precursor0.pdf','../../build/out/Average_Spline/spline/scatter/Precursor0.pdf','../../build/out/Average_Spline/spline/res/Precursor0.pdf','../../build/out/Average_Spline/spline/gof/Precursor0.txt','../../build/out/Average_Spline/spline/model/Precursor0.xml','../../build/out/Average_Spline/spline/eval/Precursor0.tab')
```
To create the model for the M+1 isotope and peptides containing any number of sulfurs:
```Matlab
IsotopeSpline(1000,'Average_spline','1','../../build/out/Average_Spline/data/Precursor1.tab','../../build/out/Average_Spline/spline/hist/Precursor1.pdf','../../build/out/Average_Spline/spline/scatter/Precursor1.pdf','../../build/out/Average_Spline/spline/res/Precursor1.pdf','../../build/out/Average_Spline/spline/gof/Precursor1.txt','../../build/out/Average_Spline/spline/model/Precursor1.xml','../../build/out/Average_Spline/spline/eval/Precursor1.tab')
```
To create the model for the monoisotope and peptides containing 0 sulfurs:
```Matlab
IsotopeSpline(1000,'0','0','../../build/out/S0/data/Precursor0.tab','../../build/out/S0/spline/hist/Precursor0.pdf','../../build/out/S0/spline/scatter/Precursor0.pdf','../../build/out/S0/spline/res/Precursor0.pdf','../../build/out/S0/spline/gof/Precursor0.txt','../../build/out/S0/spline/model/Precursor0.xml','../../build/out/S0/spline/eval/Precursor0.tab')
```
The previous commands create separate spline model .xml files. Use this command to merge them into a singe .xml file.

Usage: combineModels.py path_to_spline_xmls max_isotope_depth max_sulfur
```ShellSession
$ python ../scripts/training/combineModels.py out/ 5 5 > out/IsotopeSplines.xml
```
The publication version of IsotopeSplines.xml is already included in the IsotopeSplines branch of OpenMS that you checked out earlier. It is also in the misc/ directory.

### Figures S-1 and S-2

USAGE: plotModel.R path_to_eval_data data_dir isotope num_sulfur out_path_for_figure max_mass

Plot average spline for the monoisotope
```ShellSession
$ Rscript ../scripts/Training/plotModel.R out/spline_eval/ out/ 0 -1 out/Average_precursor0_model.eps 10000
```
Plot average spline for the M+1 isotope
```ShellSession
$ Rscript ../scripts/Training/plotModel.R out/spline_eval/ out/ 1 -1 out/Average_precursor1_model.eps 10000
```
Plot sulfur-specific spline for the monoisotope and 0 sulfurs
```ShellSession
$ Rscript ../scripts/Training/plotModel.R out/spline_eval/ out/ 0 0 out/precursor0_model.eps 10000
```

The figures are out/Average_precursor0_model.eps, etc.

### Figure S-3

USAGE: SpeedTest min_mass max_mass max_sulfurs num_samples
```ShellSession
$ ./SpeedTest 400 9500 5 100000 > out/runtimes.out
$ Rscript ../scripts/theoretical/plotRuntimeComparisons.R out/runtimes.out out/runtimes.eps
```

Figure S-3 is out/runtimes.eps

### Figure 1

USAGE: plotModelToProteome.R path_to_spline_evals out_path isotope max_sulfurs out_path_for_figure max_mass
```ShellSession
$ ./GenerateTrainingData ../data/human_sp_112816.fasta out/proteome/ 10000 5 1
$ Rscript ../scripts/theoretical/plotModelToProteome.R out/spline_eval/ out/ 0 1 out/spline_comparison_0_model.eps 10000
```

Figure 1 is out/spline_comparison_0_model.eps. Repeat for other isotopes by changing the isotope parameter.

### Figure 2 and Table S1

Here we're only going to perform 1/1000 of the simulated fragments because it takes a long time otherwise. We split it on a compute cluster for the publication.

```ShellSession
$ ./CompareToTheoretical
USAGE: CompareToTheoretical fasta_path job_id num_jobs do_frag residual_file score_file stats_file bin_size_chi bin_size_res

$ ./CompareToTheoretical ../data/human_sp_112816.fasta 1 1000 1 out/residuals_fragment.out out/scores_fragment.out out/stats_fragment.out 0.1 0.0025
$ Rscript ../scripts/theoretical/plotComparisons.R out/fragment_scores.txt out/fragment_chisquared.eps 0.1 T
```

Figure 2 is out/fragment_chisquared.eps

The file out/stats_fragment_1.txt contains the results used for Table S1.

### Figure 3

```ShellSession
$ ./CompareToTargeted ../data/Neuro_04.mzML ../data/Neuro_04_centroid.mzML out/out04.tab out/calc_out04.tab out/scores_out04.tab
$ Rscript ../scripts/experimental/lowThroughput/IndividualSpectrumIsotopes.R out/out04.tab out/calc_out04.tab out/scores_out04.tab out/low_throughput.eps
```
Figure 3 is out/low_throughput.eps

### Figure 4 and Table 2

```ShellSession
$ ./CompareToShotgun data/HELA_2017-10-25_CID25_OT.mzML data/HELA_2017-10-25_CID25_OT.idXML 0.0 out/ MS2 MS2 CID_25
$ Rscript ../scripts/experimental/highThroughput/plotShotgunResults.R out/distributionScores.out out/
```

Figure 4 is out/chi-squared_incomplete_2.pdf and out/chi-squared_incomplete_3.pdf

The values for Table 1 are printed to the console.
