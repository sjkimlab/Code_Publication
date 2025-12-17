%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last update date: 12/16/2025

Description: this script make plots of FISH data. It show how various
quantites change over different timePoints
    1) percentage of cells with signals
    2) spot number per cell
    3) spot amplitude
    4) mean xNorm value

-------- please go to the directory that contains FISH Data folder --------
                  i.e.  'yourPath\FISHanalysis2'
---------------------------------------------------------------------------
%}    

%% plotting of multiple strains
    
clear, clc, close all

fishPath = fullfile( pwd, 'FISH Data');

listStrain = dir( fullfile( fishPath, '*SK*')); % strain folders
if isempty( listStrain)
    error( '       Please change the current working path to ''FISHanalysis2'' that contains ''FISH Data folder''')
end
plotNum = getPlotNum( listStrain); % find which strains to plot 

%  'min' for SK390 & 435, SK98 & 575, 's' for SK519 & 591
sList = { 'SK519', 'SK591'}; % these strains has second as unit

signalFlag = 1; % true: cells with signal, false: include cells with no signal


for i =  plotNum % different strains
    
    dataPath = fullfile( fishPath, listStrain(i).name);
    
    % manually choose what folder to plot
    listFile = dir( fullfile( dataPath, 'FISH*')); % all folders of different days    
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
    cellStat = nan( listN, 3);      % cell spots statistics
    ampStat = nan( listN, 2);       % spot amplitude statistics
    xNormStat = nan( listN, 2);     % xNorm statistics
            
        
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
                
        % calculate mean xNorm with error by bootstrap
        funMean = @(x) mean( x, 'omitnan');
        m = bootstrp( 1000, funMean, abs( xNorm));
        xNormStat( j,:) = [ mean(m), std(m)];
    end    
        
    % Plot quantity changes over time 
    tPoints = double( erase( timeAll, {'t' 'm'}));
    
    legtxt = sprintf( '%s, %s', strain, Date);

    % percentage of cells with 1+ spots
    figure(1)
    plot( tPoints, cellStat(:,3), 'o-', 'MarkerSize', 10, 'LineWidth', 2,...
        'DisplayName', legtxt), hold on

    % number of spots/cell
    figure(2)
    errorbar( tPoints, cellStat(:,1), cellStat(:,2), 'LineWidth', 2,...
        'CapSize', 10, 'DisplayName', legtxt), hold on
    
    % signal amplitude
    figure(3)
    errorbar( tPoints, ampStat(:,1)*100, ampStat(:,2)*100, 'LineWidth', 2,...
        'CapSize', 10, 'DisplayName', legtxt), hold on
    
    % mean xNorm value
    figure(4)
    errorbar( tPoints(2:end), xNormStat(2:end,1), xNormStat(2:end,2), 'LineWidth', 2,...
        'CapSize', 10, 'DisplayName', legtxt), hold on

        
    fprintf( '\n')
end

%% Plot Setting for Time Course

% number of spots/cell
figure( 1)
set( gcf, 'Position', [200 400 350 320])
set( gca, 'FontSize', 13)
xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
ylabel( '% of cell with spots', 'FontSize', 15)
title( '% of cell with spots', 'FontSize', 15)
legend( 'Location', 'southeast', 'FontSize', 12)
xlim( [0 tmax])

% percentage of cells with 1+ spots
figure( 2)
set( gcf, 'Position', [560 400 350 320])
set( gca, 'FontSize', 13)
xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
ylabel( '# of spots/cell', 'FontSize', 15)
title( ['spot N/cell' extra], 'FontSize', 15)
legend( 'Location', 'southeast', 'FontSize', 12)
ylim( [1 2])

% amplitude 
figure(3)
set( gcf, 'Position', [920 400 350 320])
set( gca, 'FontSize', 13)
xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
ylabel( 'Spot Intensity (a.u.)', 'FontSize', 15)
title( 'Spot Intensity', 'FontSize', 15)
legend( 'Location', 'southeast', 'FontSize', 12)
xlim( [0 tmax])    

% mean xNorm value
figure(4)
set( gcf, 'Position', [1280 400 350 320])
set( gca, 'FontSize', 13)
xlabel( sprintf( 'Time (%s)', timeUnit), 'FontSize', 15)
ylabel( 'Normalized X Position', 'FontSize', 15)
title( 'xNorm', 'FontSize', 15)
legend( 'Location', 'northeast', 'FontSize', 12)
xlim( [0 tmax])
ylim( [0.25 0.5])

