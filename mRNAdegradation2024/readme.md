## About
This folder contains analysis code used in Kim et al on "Spatial and genetic constraints govern coupling between transcription, translation, and mRNA degradation in bacteria".

## FISH analysis (tested on MATLAB version R2023a)
We provide three versions of FISH analysis code depending on the need.

- `FISHanalysis2` folder contains an example dataset to test plotting. The scripts were lightly modified from those in `FISHanalysis` to use these data files accordingly.

- `FISHanalysis` folder is the original version without an example data set. This contains two additional script files that are used to 1) prepare raw files for oufti & uTrack analysis and 2) combine multiple analysis results, respectively.

- A more complete workflow for handling raw image files can be found in the [FISH2025](https://github.com/sjkimlab/FISH2025.git) repository.

## Rif-seq analysis (tested on R version 4.4.3 (R Foundation))
We provide R-code performing piece wise linear rigression to Rif-seq time course data.
The raw timepoint data is provided in an excel file. They are RNA signal (AU) at t = 0, 0.5, 1, 2, 3, 4, 5, 6, 8, 10 min after addition of Rifampicin. This data also appears in Source Data 4 file.

## SEnd-seq, Rif-seq, Ribo-seq comparison plots (tested on MATLAB version R2021a)
We provide matlab code used to plot dot density plots of PF vs TE or of half-life vs TE, as well as to calculate correlation values and Fisher test.
Source Data 3 and 4 are required to run the code.

## E. coli gene expression TASEP simulation GUI (tested on MATLAB version R2023a)
We provide matlab GUI used to run stochastic simulations and to show simulation results (data viewer).

## ODE-RK4 (tested on python version 3.8.20)
Python scripts for drawing ODE models for lacZ transient induction results (initial rise and decay) depending on kd1 and kd2 values.
The results are shown in Extended Data Figure 3a.
