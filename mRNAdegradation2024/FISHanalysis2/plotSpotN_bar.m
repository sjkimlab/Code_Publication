%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 11/28/2023
    Last update date: 12/16/2025

Description: this script is for plotting bar graph of FISH data spot number
per cell at specific time point

-------- please go to the directory that contains FISH Data folder --------
                  i.e.  'yourPath\FISHanalysis2'
---------------------------------------------------------------------------
%}    

clear, clc, close all

fishPath = fullfile( pwd, 'FISH Data');

listStrain = dir( fullfile( fishPath, '*SK*')); % strain folders
if isempty( listStrain)
    error( '       Please change the current working path to ''FISHanalysis2'' that contains ''FISH Data folder''')
end
plotNum = getPlotNum( listStrain); % find which strains to plot 

signalFlag = 0; % true: cells with signal, false: include cells with no signal

colorList = get( gca,'colororder'); close, c = 0;  

dataPath = fullfile( listStrain(1).folder, listStrain(1).name);
listFile = dir( fullfile( dataPath, 'FISH*')); % all folders of different days    
listN = length( listFile);

for t = 2: listN % timePoint

    c = c + 1;    figure( 'Position', [100+290*c 400 280 320])
    
    Nstat = cell( length( plotNum),1);  d = 0;
    
    for j = plotNum % strains
        
        dataPath = fullfile( listStrain(j).folder, listStrain(j).name);
        
        % manually choose what folder to plot
        listFile = dir( fullfile( dataPath, 'FISH*')); % all folders of different days    

        % load FISH result file of the corresponding timePoint
        load( [ listFile( t).folder '\' listFile( t).name])
        fprintf( ' ~~~ ''%s'' is Loaded ~~~\n', listFile( t).name)

        % ~~~ condition for cell variables ~~~~
        goodCells = 1: totalCells;  extra = ' (all cells)'; 
        if logical( signalFlag)
            goodCells = find( cellSpots > 0);   extra = ' (1+ spot)'; 
        end
        cSpots = cellSpots( goodCells); % number of spots (>0) in cells 
        
        % bootstrap to get errorbar
        nReps = 1000;    spotN = 0: 5;  N = nan( nReps, length( spotN));        
        for k = 1: nReps
            p = unidrnd( length( cSpots), length( cSpots), 1); % sample with replacement
            sample = cSpots( p);
            [N(k,:), ~] = histcounts( sample, 'binMethod', 'integers', 'binLimit', [-0.5 5.5], 'Normalization', 'probability');
        end
        
        % plot bar graph
        d = d + 1;  
        bar( spotN, mean(N), 1, 'FaceColor', colorList( d,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w', 'LineWidth', 1.5, ...
            'DisplayName', strain), hold on % sprintf( '%s, %d cells', strain, length( cSpots))
                
        errorbar( spotN, mean( N), std( N), 'Color', [0 0 0], 'CapSize', 10,...
            'LineWidth', 1, 'LineStyle', 'none', 'HandleVisibility', 'off') 
        
        Nstat{ d} = [ mean(N); std(N)]; % statistics of the bar percentage        
    end
    
    
    % Plot Setting
    hold off
    figure( gcf)
    set( gca, 'FontSize', 14)
    xlabel( 'Number of spots per cell', 'FontSize', 15)
    ylabel( 'Probability', 'FontSize', 15)
    title( sprintf( '%s, %s', Date, timePoint), 'FontSize', 15)
    legend( 'Location', 'northeast', 'FontSize', 13)

    if ~signalFlag
        xlim( [-0.5 4.8])
    else
        xlim( [0.5 4.8])
    end
    ylim( [0 0.6])
        
    fprintf( '\n')
end

