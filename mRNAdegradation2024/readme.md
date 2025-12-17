## About
This folder contains analysis code used in Kim et al on mRNA degradation in bacteria.

## FISH analysis
We provide three versions of FISH analysis code, depending on the need.

- `FISHanalysis2` folder contains an example dataset to test plotting. The scripts were lightly modified to use these data files accordingly.

- `FISHanalysis` folder is the original version without an example data set. This contains two additional script files that are used to 1) prepare raw files for oufti & uTrack analysis and 2) combine multiple analysis results, respectively.

- A more complete workflow for handling raw image files can be found in the [FISH2025](https://github.com/sjkimlab/FISH2025.git) repository.

  

## ODE-RK4
Python scripts for drawing ODE models for lacZ transient induction results (initial rise and decay) depending on kd1 and kd2 values.
The results are shown in Extended Data Figure 3a.

