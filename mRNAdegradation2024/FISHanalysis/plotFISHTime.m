%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last updated at 4/2/2024

Description: this script make plots of FISH data. It show how various
quantites change over different timePoints
    1) percentage of cells with signals
    2) spot number per cell
    3) spot amplitude
    4) mean xNorm value (optional)

---------------------------------------------------------------------------
%}    

%% plotting of multiple strains
    
clear, clc, close all

fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\FISH\';

listStrain = dir( [fishPath '*SK*']); % strain folders
plotNum = getPlotNum( listStrain); % find which strains to plot 

%  'min' for SK390 & 435, SK98 & 575, 's' for SK519 & 591
sList = { 'SK519', 'SK591'}; % these strains has second as unit

txtToInclude = 'comb123'; % what folder do you want to plot

signalFlag = 1; % true: cells with signal, false: include cells with no signal
xNormMeanFlag = 0; % flag: plot mean xNorm over time


for i =  plotNum % different strains
    
    dataPath = [ fishPath listStrain(i).name '\'];
    
    % manually choose what folder to plot (optional)
%     list = dir( [dataPath '2*']); % all folders of different days
%     listDate = list( [list.isdir]); % only look at folders instead of files        
%     plotDay = getPlotNum( listDate); % find which day of data you want to plot
%     list = dir( [dataPath listDate( plotDay).name '\FISH*']);

    list = dir( [dataPath '*' txtToInclude '*\FISH*']); % find files to plot
    listFile = sortNat( list); % sort timePoints by natural order    
    
    % ~~~ Plot only part of the data ~~~~
%     listFile = listFile( 2:4); % timePoints to plot
    listN = length( listFile); % number of timePoints
    
    load( [ listFile( 1).folder '\' listFile( 1).name])    
    fprintf( ' ~~~~~~~~~  ''%s-%s''  Loaded ~~~~~~~~~\n', strain, Date)
    
    if ismember( strain, sList) 
        % find the right unit for time
        timeUnit = 's';     tmax = 170; % max value for x-axis in plotting
    else
        timeUnit = 'min';   tmax = 6; 
    end
    
    % initialize variables for multiple T
    timeAll = strings( listN, 1);   % stores timePoints
    cellStat = nan( listN, 3);  % cell spots statistics
    ampStat = nan( listN, 2);   % spot amplitude statistics
    xNormStat = nan( listN, 2); % xNorm statistics
            
        
    for j = 1: listN % timePoints
        
        load( [ listFile( j).folder '\' listFile( j).name]) % load file
        
        timeAll( j) = timePoint; % string
        
        % ~~~ condition for cell variables ~~~~
        goodCells = 1: totalCells;  extra = ' (all cells)'; 
        if logical( signalFlag)
            goodCells = find( cellSpots > 0);   extra = ' (1+ spot)'; 
        end
        
        condSpots = ismember( spotCell, goodCells); % flag for spots in selected cells
        xNorm = spotNorm( condSpots, 1);
        amp = spotAmp( condSpots);
        ampStat( j,:) = [ mean( amp), std( amp)/ sqrt( length(amp))];
                
        % calculate cell statistics (how many spots, % of cell with spots)
        cSpots = cellSpots( goodCells); cellPerc = 100* sum( cellSpots>0)/ totalCells;
        cellStat( j,:) = [mean( cSpots), std( cSpots)/ sqrt( length(cSpots)), cellPerc];
                        
        % display # of cells & spots used in analysis
        fprintf( ' ~~~ %s %4s: %6d spots, %6d cells has signal ~~~\n',...
            strain, timePoint, sum( cellSpots), sum( cellSpots>0))
                
        if logical( xNormMeanFlag)
            % calculate mean xNorm with error by bootstrap
            funMean = @(x) mean( x, 'omitnan');
            m = bootstrp( 1000, funMean, abs( xNorm));
            xNormStat( j,:) = [ mean(m), std(m)];
        end        
    end    
        
    % Plot quantity changes over time 
    tPoints = double( erase( timeAll, 't'));
%     Date = 'text'; % for display
    
    % percentage of cells with 1+ spots
    figure(10)
    plot( tPoints, cellStat(:,3), 'o-', 'MarkerSize', 10, 'LineWidth', 2,...
        'DisplayName', sprintf( '%s, %s', strain, Date)), hold on

    % number of spots/cell
    figure(11)
    errorbar( tPoints, cellStat(:,1), cellStat(:,2), 'LineWidth', 2,...
        'CapSize', 10, 'DisplayName', sprintf( '%s, %s', strain, Date)), hold on
    
    % signal amplitude
    figure(12)
    errorbar( tPoints, ampStat(:,1)*100, ampStat(:,2)*100, 'LineWidth', 2,...
        'CapSize', 10, 'DisplayName', sprintf( '%s, %s', strain, Date)), hold on
    
    if logical( xNormMeanFlag)
        figure(13)
        errorbar( tPoints, xNormStat(:,1), xNormStat(:,2), 'LineWidth', 2,...
            'CapSize', 10, 'DisplayName', sprintf( '%s', strain)), hold on
    end
        
    fprintf( '\n')
end

%% Plot Setting for Time Course

% number of spots/cell
figure( 10)
set( gcf, 'Position', [200 400 350 320])
set( gca, 'FontSize', 13)
xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
ylabel( '% of cell with spots', 'FontSize', 15)
title( '% of cell with spots', 'FontSize', 15)
legend( 'Location', 'southeast', 'FontSize', 13)
xlim( [0 tmax])

% percentage of cells with 1+ spots
figure( 11)
set( gcf, 'Position', [560 400 350 320])
set( gca, 'FontSize', 13)
xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
ylabel( '# of spots/cell', 'FontSize', 15)
title( ['spot N/cell' extra], 'FontSize', 15)
legend( 'Location', 'southeast', 'FontSize', 13)
ylim( [1 2])

% amplitude 
figure(12)
set( gcf, 'Position', [920 400 350 320])
set( gca, 'FontSize', 13)
xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
ylabel( 'Spot Intensity (a.u.)', 'FontSize', 15)
title( 'Spot Intensity', 'FontSize', 16)
legend( 'Location', 'southeast', 'FontSize', 13)
xlim( [0 tmax])    

if logical( xNormMeanFlag)
    figure(13)
    set( gcf, 'Position', [1280 400 350 320])
    set( gca, 'FontSize', 13)
    xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
    ylabel( 'Normalized X Position', 'FontSize', 15)
    title( 'xNorm', 'FontSize', 16)
    legend( 'Location', 'northeast', 'FontSize', 13)
    xlim( [0 tmax])
end

