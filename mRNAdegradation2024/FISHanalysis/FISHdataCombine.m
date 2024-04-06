%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 9/12/2023
    Last updated at 3/19/2024

Description: this script is for combining multiple days FISH data

---------------------------------------------------------------------------
%}    

%% Choose what days & timePoints of data to combine

clear, clc, close all

combDate = '240403 comb';

fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\FISH\';

strainToComb = input( 'What strain do you want to combine (like 519): ', 's');

dataPath = [ fishPath 'SK' strainToComb '\'];


% Determine what days of data to combine
list = dir( [dataPath '2*']);
listDate = list( [list.isdir]); % only look at folders instead of files

    fprintf( '\n~~~~~ There are %d files for this strain ~~~~~\n', length( listDate));
    for k = 1: length( listDate)
        disp( ['  ' num2str(k) '. ' listDate(k).name])
    end
    tmp = string( input( '\nWhich days of data do you want to combine(like: 1 3 4): ', 's')); % just input: 1 3 4 5
    combDay = double( regexp( tmp, '\d+', 'match')); % in case you input extra space
        
% Determine what timePoints of data to combine
list = dir( [dataPath listDate( combDay(1)).name '\FISH*']); % list FISH files of all timePoints
listFile = sortNat( list);

    fprintf( '\n~~~~~ There are %d files for this strain ~~~~~\n', length( listFile));
    for k = 1: length( listFile)
        disp( ['  ' num2str(k) '. ' listFile( k).name])
    end
%     tmp = string( input( '\nWhich files (timePoints) do you want to combine (like: 1 3 4): ', 's'));
%     combNum = double( regexp( tmp, '\d+', 'match')); % in case you input extra space
    combNum = 1: length( listFile);
        
%% Combine files from different days

for i = combNum 
    
    spotCellComb = [];  spotNormComb = [];  spotAmpComb = []; 
    cellAreaComb = [];  cellLengthComb = [];    cellSpotsComb = [];  cellSpotStatComb = {};

    for j = combDay
        
        listFile = sortNat( dir( [dataPath listDate( j).name '\FISH*']));
        load( [listFile( i).folder '\' listFile( i).name])
        fprintf( '\n  ~~~~ %s Loaded  ~~~~', listFile( i).name)
                
        % store the information from this image and stack them together,
        spotCellComb = [ spotCellComb; spotCell]; % cell number assigned to spots (accumulated)
        spotNormComb = [ spotNormComb; spotNorm]; % spot normalized position
        spotAmpComb = [ spotAmpComb; spotAmp];  % added at 9/20/2023

        cellAreaComb = [ cellAreaComb; cellArea   ];
        cellLengthComb = [ cellLengthComb; cellLength];
        cellSpotsComb = [ cellSpotsComb; cellSpots]; % # of spots in all cells
        
        % record the origin of the files
        cellSpotStat(:,2) = { Date};
        cellSpotStat(:,3) = { timePoint};
        cellSpotStat(:,4) = { listFile( i)};
        cellSpotStatComb = [ cellSpotStatComb; cellSpotStat]; % count each cell's spot #
    end
    
    spotCell = spotCellComb;
    spotNorm = spotNormComb;
    spotAmp = spotAmpComb;
    
    cellArea = cellAreaComb;
    cellLength = cellLengthComb;
    cellSpots = cellSpotsComb;
    cellSpotStat = cellSpotStatComb;
    
    combRecord = unique( cellSpotStat(:,2));
    
    cellWid = cellArea./ cellLength;  % unit: um
    totalCells = length( cellWid);
    totalSpots = length( spotAmp);
    
    savePath = [ fishPath strain '\' combDate '\'];
    
    if ~exist( savePath, 'dir')
        mkdir( savePath)
    end
    
    Date = combDate;
    fishName = [ 'FISH ' strain ' ' timePoint ' ' Date];
    save( [ savePath fishName], 'pixelSize', 'totalSpots', 'totalCells', 'combRecord', ...
        'cellArea', 'cellLength', 'cellWid', 'cellSpots', 'cellSpotStat', ...
        'spotCell', 'spotNorm', 'spotAmp', 'Date', 'strain', 'timePoint', 'savePath', 'fishName')    
    
    fprintf( '\n~~~~~~ %s Saved ~~~~~~\n', [fishName '.mat'])
end
    
% cd( savePath)
