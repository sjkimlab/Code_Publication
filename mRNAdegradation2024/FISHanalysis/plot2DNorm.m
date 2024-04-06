%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last updated at 4/2/2024

Description: this script make plots for FISH data. It plots the 2D Norm
histogram over different timePoint. 

FISHdataAnalysis output:
    spotCell: cell number for each spots if it's inside any cell
    spotNorm: [xNorm, lNorm] within the cell for each spots
    cellSpots: the number of spots in each cells, 0 if no signal
    cellArea & cellLength & cellWid: properties of each cells (unit: um)

---------------------------------------------------------------------------
%}    

%% plotting of multiple strains
    
clear, clc, close all

fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\FISH\';

listStrain = dir( [fishPath '*SK*']); % strain folders
plotNum = getPlotNum( listStrain); % find which strains to plot 

%  'min' for SK390 & 435, SK98 & 575, 's' for SK519 & 591
sList = { 'SK519', 'SK591'}; % these strains has second as unit

% txtToInclude = 'comb123'; % what folder do you want to plot

signalFlag = 1; % true: cells with signal, false: include cells with no signal
c = 0;

% parameters for plotting
norm2DScale = 0.05; % colormap limit
nBinY = 5;  nBinX = 3* nBinY; % number of bins for 2D Norm heatmap
maxN( 1, length( plotNum)) = 0; % stores the max value for colorbar plot


for i =  plotNum % different strains
    
    dataPath = [ fishPath listStrain(i).name '\'];
    
    % manually choose what folder to plot (optional)
    list = dir( [dataPath '2*']); % all folders of different days
    listDate = list( [list.isdir]); % only look at folders instead of files        
    plotDay = getPlotNum( listDate); % find which day of data you want to plot
    list = dir( [dataPath listDate( plotDay).name '\FISH*']);

    list = dir( [dataPath '*' txtToInclude '*\FISH*']); % find files to plot
    listFile = sortNat( list); % sort timePoints by natural order    
    
    % ~~~ Plot only part of the data ~~~~
    listFile = listFile( 2:4); % timePoints to plot
    listN = length( listFile); % number of timePoints
    
    load( [ listFile( 1).folder '\' listFile( 1).name])    
    fprintf( ' ~~~~~~~~~  ''%s-%s''  Loaded ~~~~~~~~~\n', strain, Date)
    
                
    % set up tiled layouts (subplots) for 2D Norm plot
    c = c + 1;
    figure( 'Position', [370*c 600-listN*90 360 70+140*listN])
    t = tiledlayout( listN, 1, 'TileSpacing', 'Compact', 'Padding', 'compact'); % one column
    tileN = 1: listN;
    
    
    for j = 1: listN % timePoints
        
        load( [ listFile( j).folder '\' listFile( j).name]) % load file
                        
%         % display the cell length & wid
%         fprintf( ' %s %s  Cell Length: %.2f,  wid: %.2f um\n',...
%             strain, timePoint, mean( [ cellLength( cellSpots > 0), cellWid( cellSpots > 0)]))
        
%         % find short 25% cells
%         [cLengthSort, order] = sort( cellLength);
%         smallCells = order( 1: round( end* 0.25)); % short 25% of cell length
%         goodCells = smallCells;     extra = ' (25% length)';
        
        % ~~~ condition for cell variables ~~~~
        goodCells = 1: totalCells;  extra = ' (all cells)'; 
        if logical( signalFlag)
            goodCells = find( cellSpots > 0);   extra = ' (1+ spot)'; 
        end        
        condSpots = ismember( spotCell, goodCells); % flag for spots in selected cells
        xNorm = spotNorm( condSpots, 1);
        lNorm = spotNorm( condSpots, 2);
                
        % display # of cells & spots used in analysis
        fprintf( ' ~~~ %s %4s: %6d spots, %6d cells has signal ~~~\n',...
            strain, timePoint, sum( cellSpots), sum( cellSpots>0))        
        
        
        % hist for a quarter of cell, flip twice to form whole cell
        [N, xEdges, yEdges] = histcounts2( abs( lNorm-0.5)*2, abs( xNorm)/2, [ nBinX, nBinY],...
            'XBinLimits', [0 1], 'YBinLimits', [0 0.5], 'Normalization', 'probability');
        % N: 15 rows * 5 columns, it should be flipped (transpose') when plotted by pcolor
        % N' is the bottom right part of the matrix, should flip to fill
        
        % flip: upside down,    flip(N,2): left-right
        allN = [ flip( flip(N',2))   flip(N');   
                       flip(N',2)         N' ];
        
        allN( end+1, end+1) = 0; % the color of each grid was determined by the first [X,Y] 

        [l, r] = size( N);
        X = -l: 1: l;   Y = -r: 1: r;

        nexttile( tileN(j))
        s = pcolor( X, Y, allN); hold on            

        s.LineWidth = 1; % edge line width
        set( gca, 'xtick',[], 'ytick',[]) % no tick or labels
        title( sprintf( '%s', timePoint), 'FontSize', 14)
        caxis( [0 norm2DScale]) % set colormap limit
        axis image, colormap jet

        maxN( j, c) = max( N(:));
        
        % plot cell mesh
        rMesh = r*0.95;     th = (0: pi/100 : pi) + pi/2;
        x = [ rMesh*cos( th)+ r-l,  rMesh*cos( th+pi)+ l-r r-l];
        y = [ rMesh* sin( th),  rMesh* sin( th+pi) rMesh];
        plot( x, y, 'w', 'lineWidth', 3)       
        
    end    
    
    % Plot Setting        
    title( t, sprintf( '%s %s', strain, Date), 'FontSize', 16)
    fprintf( '\n~~~  maxN for colorbar should be %.3f ~~~\n', max( maxN(:)))
    
    fprintf( '\n')
end
