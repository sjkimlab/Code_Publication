%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 2023
    Last updated at 3/19/2024

Description: this script ask you to choose files from a list for plotting

---------------------------------------------------------------------------
%}

function plotNum = getPlotNum( list)

    plotNum = 1; % for multiple file of the same strain    
    
    if length( list) > 1
        fprintf( '\n~~~~~ There are %d files for this strain ~~~~~\n', length(list));
        for k = 1:length(list)
            disp( ['  ' num2str(k) '. ' list(k).name])
        end
        tmp = string( input( '\nWhich files do you want to plot (like: 1 3 4): ', 's')); % just input: 1 3 4 5
        plotNum = double( regexp( tmp, '\d+', 'match')); % find all numbers but not extra space
    end
    
fprintf( '\n')
end

