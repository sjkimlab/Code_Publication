## About

This repo contains analysis and plotting scripts for smFISH data. An example dataset is also provided, which you can use to generate plots. The dataset already integrates outputs from **Oufti** & **u-track** analyses, and is ready for plotting. 

To run the plotting code:

1. Add the `FISHanalysis2` folder and all its subfolder to the MATLAB search path. 
2. Set the MATLAB current working directory to the `FISHanalysis2` folder.



#### System Requirement

- MATLAB (R2023a or later)

#### Analysis Code

- `FISHdataAnalysis.m`: this script combines the output from oufti & u-track analyses and performs further analysis, such as assigning foci to cells and calculating normalized cellular position.

#### Plotting Code

- `plot2DNorm.m`: plot 2D histogram of normalized spot localizations.
- `plotFISHTime.m`: plot how various quantities change over time  
  - *percentage of cells with signals, spot number per cell, spot intensity, mean xNorm*

- `plotSpotN_bar.m`: plot bar graph of spot numbers per cell.

#### Dataset

Example dataset of two strains (SK390 & SK435) is provided in the `FISH Data` folder. Each strain subfolder contains FISH data files at 6 different time points.   



## Running Workflow

1. Download MATLAB and add the `FISHanalysis2` folder and all its subfolder to MATLAB search path.
2. Set the MATLAB current working directory to the `FISHanalysis2` folder. For example, enter `cd( 'yourPath\FISHanalysis2')` in the MATLAB command window
3. Run plotting scripts to visualize analysis results
