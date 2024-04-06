%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/28/2023
    Last updated at 3/13/2024

Description: this script organized the FISH tiff files output into
corresponding folder for phase, Cy5, Cy3. Then it creates folders for Cy3
specifically with renamed files for uTrack analysis


For FISH experiments, the images taken at different time points are stored
in folders named as t0m, t1m, t2m... (for minutes) or t0, t30, t60, t90...
(for seconds). The .nd files are converted into tiff files using Nikon NIS
Elements software.

The tiff files are named with their channel numbers:
    c1 - Cy5 | c2 - Cy3 | c3 - phase
This script moves the corresponding tiff files into the right folders
    1ph | 2Cy5 | 3Cy3
---------------------------------------------------------------------------
%}


% please go to the directory that stores folders of different timePoints
%     e.g. t0, t30, t60 ... or t0m, t1m, t2m ...

%% 1. Organize tiff files in different channels

clear, clc

list = dir( 't*'); % FISH folders were named as t1m or t30
list = list( [list.isdir]); % only look at folders instead of files
folderName = { '2Cy5' '3Cy3' '1ph'}; % corresponds to channel 1, 2, 3


for i = 1: length( list)
    
    path = [ list(i).folder '\' list(i).name '\'];  % SK519\t20 or SK98\t0m
    cd( path)
    
    fprintf( '~~~    timePoint: %s    ~~~\n', list(i).name)
    
    for c = 1: 3
        
        % Split images from different channels into folders
        listTif = dir( sprintf( 'epi*c%d*', c)); % find images of this channel
        mkdir( folderName{c})
        
        for k = 1: length( listTif)
            
            name = listTif(k).name; % original name from FISH imaging (e.g. epi170c1.tif)
            movefile( name, [ folderName{c} '\' name]) % move into 2Cy5 Folder
        end
        fprintf( '~~~  %s tiff files Moved  ~~~\n\n', folderName{c})
    end
end
cd ..\


%% 2. Rename the epi Cy3 images for uTrack analysis

for i = 1: length( list)
    
    path = [ list(i).folder '\' list(i).name '\'];  % SK519\t30 or SK98\t1m
    folderName = ['Cy3_' list(i).name];             % Cy3_t30
    
    % copy the epi Cy3 channel to another folder named Cy3_time
    copyfile( [path '3Cy3'], [path '..\' folderName]);  cd( [path '..\' folderName])
    
    % rename epi files so that uTrack can treat them as a time-stack (I
    % only use this to run uTrack spotDetection more efficiently, but not
    % link them to form tracks)
    %       epi001c2.tif --> epiCy3t01.tif
    %       epi002c2.tif --> epiCy3t02.tif    
    %
    listEpi = dir( 'epi*');
    for k = 1: length( listEpi)        
        name = listEpi(k).name; % original name of FISH file (e.g. epi001c2.tif)
        newName = sprintf( 'epiCy3t%02.f.tif', k); 
        movefile( name, newName)
        fprintf( '~~~  %s -> %s  ~~~\n', name, newName)
    end
    
    fprintf( '~~~~~~    %s File Conversion Complete    ~~~~~~\n\n', list(i).name)
end
cd ..\


%% 3. Run uTrack on the FISH Cy3 data for spot detection

fprintf( '~~~ Run uTrack, select input channel folder ''Cy3_time'', set the output folder the same ~~~')
% uTrack output will be saved to the same folder containing Cy3 tiff files

% for FISH data, I'm currently using following parameter fow spot detection
%   Gaussian std = 1.28 pixel (Hamamatsu camera)
%   Alpha (Detection) = 0.015
%   Alpha (fitting) = 0.05 (default)

