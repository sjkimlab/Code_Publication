# Ecoli_geneExpression_TASEP

Ecoli_geneExpression_TASEP is a matlab application to simulate transcription, translation, and RNA degradation and to analyze the simulated datasets. Importantly, premature transcription termination is allowed when RNA polymerase and ribosome are not coupled.

This program requires either MATLAB or MATLAB Runtime, which does not require the user to have a MATLAB license.


## Using MATLAB runtime

The .exe version of the code can be run via MATLAB Runtime. You can download and install this runtime from:

`https://www.mathworks.com/products/compiler/mcr/index.html`

Our code is tested on MATLAB runtime R2023a(9.14).

After downloading the Runtime, unzip the file and double-click the setup file.


## GUI, using MATLAB App designer

The .mlapp version of the code can be run via MATLAB app designer. The .mlapp file can be found in the following repository:
`\uncompiledCode\Ecoli_geneExpression_TASEP.mlapp`

The .mlapp requires the following file to be in the same directory:
`Ecoli_geneExpression_TASEP_function.m`

We recommend setting the MATLAB Current Folder to `\uncompiledCode\` then double-clicking `Ecoli_geneExpression_TASEP_dataViewer.mlapp` in the MATLAB Current Folder window.
The App Designer start page will open, which may take some time. Press the Run button in the App Designer window to launch the GUI.

This method was tested on MATLAB R2021a and R2023a.

## Running the GUI
A more detailed guide to help use the simulation GUI can be found under the name: **`GUI_user_guide.pdf`**

