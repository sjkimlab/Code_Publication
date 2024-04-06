%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 11/28/2023
    Last updated at 3/22/2024

Description: this script is for plotting bar graph of FISH data spot number
per cell at specific time point

12.12: added bootstrap for getting a error bar for the bar graph
3.21: added option to plot gene loci on top
---------------------------------------------------------------------------
%}    
    
clear, clc, close all

fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\FISH\';
figPath = [ fishPath 'Figure\'];

listStrain = dir( [fishPath '*SK*']); %  combined strains
plotNum = getPlotNum( listStrain); % find which strains you want to plot 
      
lociFLag = 1; % flag for plotting loci number
signalFlag = 1; % true: cells with signal, false: include cells with no signal
saveFlag = 0; % flag for saving plots into files

txtToInclude = 'comb123'; % what folder do you want to plot

colorList = get( gca,'colororder'); close, c = 0;  

for t = 3: 4 % 1: 5 % timePoint

    c = c + 1;    figure( 'Position', [100+270*c 400 260 340])
    
    Nstat = cell( length( plotNum),1);  d = 0;
    
    for j = plotNum % strains
        
        dataPath = [ listStrain(j).folder '\' listStrain(j).name '\' ];
        
%         % manually choose what folder to plot (optional)
%         list = dir( [dataPath '2*']); % all folders of different days
%         listDate = list( [list.isdir]); % only look at folders instead of files        
%         plotDay = getPlotNum( listDate); % find which day of data you want to plot
%         list = dir( [dataPath listDate( plotDay).name '\FISH*']);
        
        list = dir( [dataPath '*' txtToInclude '*\FISH*']);
        listFile = sortNat( list); % sort timePoints by natural order

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
        bar( spotN, mean(N), 0.9, 'FaceColor', colorList( d,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w', 'LineWidth', 1.5, ...
            'DisplayName', strain), hold on % sprintf( '%s, %d cells', strain, length( cSpots))
                
        errorbar( spotN, mean( N), std( N), 'Color', [0 0 0], 'CapSize', 10,...
            'LineWidth', 1, 'LineStyle', 'none', 'HandleVisibility', 'off') 
        
        Nstat{ d} = [ mean(N); std(N)]; % statistics of the bar percentage        
    end
    
    if lociFLag 
        % plot SK71 gene loci number
        load( 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\2024\240127-SK71 SJK loci\spotNperCell.mat', 'cellSpots') % where I store my loci data
        cSpots = cellSpots( cellSpots > 0); % ignore cells with no loci signal
        [Nloci, ~] = histcounts( cSpots, 'binMethod', 'integers', 'binLimit', [-0.5 5.5], 'Normalization', 'probability');
        scatter( spotN, Nloci, 60, 'r', 'LineWidth', 2.5, 'DisplayName', 'loci')
    end
    
    % Plot Setting
    hold off
    figure( gcf)
    set( gca, 'FontSize', 13)
    xlabel( 'Number of spots per cell', 'FontSize', 15)
    ylabel( 'Probability', 'FontSize', 15)
    title( sprintf( '%s, %s', Date, timePoint), 'FontSize', 15)
    legend( 'Location', 'northeast', 'FontSize', 13)

    if ~signalFlag
        xlim( [-0.5 4.8])
    else
        xlim( [0.5 4.8])
    end
    ylim( [0 0.65])
        

    if logical( saveFlag)
%         savefig( [figPath sprintf( 'spotNumber SK390 & 435 at %s', timePoint)])
%         savefig( [figPath sprintf( 'spotNumber SK575 & 98 at %s', timePoint)])
        cd( figPath)
    end
    
    fprintf( '\n')
end

