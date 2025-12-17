%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last updated at 3/19/2024

Description: this script is for FISH data analysis, please run it section
by section 

====== Output =======
spotCell: cell number for each spots if it's inside any cell, nan for none
spotNorm: [xNorm, lNorm] within the cell for each spots, [nan nan] for none
spotAmp: signal amplitude of all spots
cellSpots: the number of spots in each cells, 0 if no signal
cellSpotsStat: unique cell ID vs number of spots
cellArea & cellLength & cellWid: properties of each cells (unit: um)

8/27/2023: added the part for dealing with no cell in a image
9/19/2023: added spot intensity as spotAmp
12/18/2023: deal with slopes with opposite sign issue
3/14/2024: moved uTrack preparation part to other script
---------------------------------------------------------------------------
%}


%% 1. Combine uTrack & oufti output and save into files 

clear, clc

% Date = '240402';    strain = 'SK591';

% specify the strain name  
Date = input( ' Please input the Date of the experiment (like: 230601):  ', 's'); fprintf( '\n')
strain = input( ' Please input the strain number (like: SK98):  ', 's'); fprintf( '\n')

% specify the folder location of FISH analysis files
fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\FISH\';
savePath = [ fishPath strain '\' Date '\'];

if ~exist( savePath, 'dir')
    mkdir( savePath)
end

listFolder = sortNat( dir( 'Cy3_t*')); % uTrack output folders

for i = 1: length( listFolder)
    
    % load uTrack output - movieInfo (containes spot coordinates)
    path = [ listFolder(i).folder '\' listFolder(i).name '\'];
    load( [path '\TrackingPackage\GaussianMixtureModels\detections_for_channel_1\Channel_1_detection_result'],...
        'movieInfo')
    
    tmp = split( listFolder(i).name, '_');  % e.g. 'Cy3_t30'
    timePoint = tmp{2};                     % e.g. 't30'
    
    % load oufti mesh output (cell meshes, e.g. 't30-mesha.mat')
    load( [ timePoint '-mesha.mat'], 'cellList', 'cellListN')
    
    % set up cell mesh and remove empty cells
    cellMeshAll = setCellMesh( cellList, cellListN);
    cellListN = cellfun( @length, cellMeshAll); % number of cell in each image
    
    totalSpots = sum( cellfun( @length, {movieInfo.amp}'));
    totalCells = sum( cellListN);
    fprintf( '~~~  %8s | spot number: %5d,  cell number: %5d', listFolder(i).name, totalSpots, totalCells)
    
    % save uTrack spot output 'movieInfo' & other info into Cy3_spotsMesh_time
    save( [ 'Cy3_spotsMesh_' timePoint], '-regexp', '^(?!(i|tmp|cellList)$).') % save variables except
    
    fprintf( '   ~~ spotsMesh Saved ~~~\n')
end

fprintf( '\n~~~~~~ uTrack & oufti output extracted and saved to Cy3_spotsMesh ~~~~~~\n\n')


%% 2. Calculate the spotNorm & other quantities, Save to files

clear

% find folder order by timePoints
list = sortNat( dir( 'Cy3_spotsMesh_*.mat'));

pixelSize = 64.5; % unit: nm

for k = 1: length( list) % timePoints
    
    load( list( k).name)
        
    spotCellAll = {};   spotNormAll = {};   spotAmpAll = {};
    cellArea = [];      cellLength = [];    cellSpots = [];   cellSpotStat = [];
    
    cellNcum = cumsum( [ 0; cellListN]); % this is for assigning unique cell number

    for i = 1: length( movieInfo) % image number
        
        cellMesh = cellMeshAll{ i};
        if isempty( fieldnames( cellMesh))
            continue % skip this image if no cells are detected
        end
        
        spotPos = [movieInfo(i).xCoord(:,1) movieInfo(i).yCoord(:,1)];
        spotAmp = movieInfo(i).amp(:,1); % added @9/19/2023
                
        nSpots = length( spotAmp);  % total number of spots detected
        spotCell = nan( nSpots, 1); % which cell the spots are in
        spotNorm = nan( nSpots, 2); % normalized position of spots [lNorm, xNorm]
        
        for Cell = 1: cellListN( i) % number of cells in this image            
            
            % a. find which cell the spots belong to
            meshOut = double( cellMesh( Cell).meshOut);
                % compare spot positions & mesh to find which cells they belong to
                inCell = inpolygon( spotPos(:,1), spotPos(:,2), meshOut(:,1), meshOut(:,2)); % return 0 or 1 in a vector
                spotCell( inCell) = Cell; % assign the spots to cell
            
            % b. find normalized position of the spots in cells
            badCell = false;
            for spotNum = find( inCell)'
                pt = spotPos( spotNum, :); % spot coordinate
                [ spotNorm( spotNum, :), badCell] = findNormPos( pt, cellMesh( Cell), false);
            end
            
            if badCell                
                warning( '   ~~~ Movie #%d Cell #%d mesh curvature has problem ~~~\n', i, cellMesh( Cell).cellId)
            end
            if cellMesh( Cell).area/ cellMesh( Cell).length * pixelSize/1000 < 0.3
                warning( '   ~~~ Movie #%d Cell #%d width is abnormal (too thin) ~~~\n', i, cellMesh( Cell).cellId)
            end
        end
        
        % count the # of spots in each cell (integer bins)
        [N, ~] = histcounts( spotCell, 'binLimit', [0.5 length( cellMesh)+0.5],...
            'binMethod', 'integers');
    %         histogram( N)

        % store the information from this image and stack them together,
        % cell number adds up over images
        spotCellAll = [ spotCellAll; spotCell + cellNcum(i)]; % cell number assigned to spots (accumulated)
        spotNormAll = [ spotNormAll; spotNorm]; % spot normalized position
        spotAmpAll = [ spotAmpAll; spotAmp]; % spot intensity
        cellArea = [ cellArea; [cellMesh.area]'];
        cellLength = [cellLength; [cellMesh.length]'];
        cellSpots = [ cellSpots; N'];        % # of spots in all cells
        cellSpotStat{i,1} = [ [cellMesh.cellId]' N']; % count each cell's spot number for all movies
    end

    cellArea = cellArea* pixelSize^2/ 10^6; % unit: um^2
    cellLength = cellLength * pixelSize/ 1000; % unit: um
    cellWid = cellArea./ cellLength;  % unit: um
    
    spotCell = cell2mat( spotCellAll);
    spotNorm = cell2mat( spotNormAll);
    spotAmp  = cell2mat( spotAmpAll);
    
%     save( ['spotResult_' timePoint], 'pixelSize', 'totalSpots', 'totalCells', ... % % 'cellMeshAll', 'cellListN', 
%         'cellArea', 'cellLength', 'cellWid', 'cellSpots', 'cellSpotStat', ... % 'spotCellAll', 'spotNormAll', 'spotAmpAll',
%         'spotCell', 'spotNorm', 'spotAmp', 'Date', 'strain', 'timePoint')
    
    fishName = [ 'FISH ' strain ' ' timePoint ' ' Date];
    save( [savePath fishName],  'pixelSize', 'totalSpots', 'totalCells', ...
        'cellArea', 'cellLength', 'cellWid', 'cellSpots', 'cellSpotStat', ...
        'spotCell', 'spotNorm', 'spotAmp', 'Date', 'strain', 'timePoint', 'savePath', 'fishName')
    fprintf( '~~~  timePoint %4s calculation completed & saved  ~~~\n', timePoint)
end

fprintf( '\n~~~~~~ Analysis Done ~~~~~~\n\n')


%% Function

function cellMeshAll = setCellMesh( cellList, cellListN)
% this function set up cell mesh and remove empty cells

    cellMeshAll = {};
    for round = 1: size( cellListN, 2) % number of images

        clear cellMesh

        nCells = size( cellList.meshData{ round}, 2);
        if nCells == 0  % in case no cell are detected
            cellMesh( 1, 1) = struct;
            cellMeshAll = [cellMeshAll; cellMesh];
            continue
        end
        cellMesh( nCells, 1) = struct;  badCells = [];
        
        for Cell = 1: nCells
            
            mesh = cellList.meshData{ round}{ Cell}.mesh; % cell outline: [n,4] - (x1, y1, x2, y2), (right, left)

            if length( mesh) > 10 % some cellMesh only have 6 points (maybe error?)
                meshMid = [ mean( mesh( :, [1 3]), 2), mean( mesh( :, [2 4]), 2)]; % midline along the long axis
                % save the cell mesh info for later use
                cellMesh( Cell).mesh = mesh;
                cellMesh( Cell).meshOut = [ mesh(:, 1:2); flipud( mesh(:, 3:4))]; % reshape the mesh matrix to form a circle [2n, 2]
                cellMesh( Cell).meshMid = meshMid;
                cellMesh( Cell).gridLen = vecnorm( diff( meshMid), 2, 2); % length of each grid (L direction)
                cellMesh( Cell).gridLenCum = cumsum( vecnorm( diff( meshMid), 2, 2)); % cumulative length of each grid (L direction)
                cellMesh( Cell).area = double( polyarea( cellMesh( Cell).meshOut(:,1), cellMesh( Cell).meshOut(:,2)));
                cellMesh( Cell).length = double( sum( cellMesh( Cell).gridLen));
                if Cell > length( cellList.cellId{ round}) % sometimes oufti has this error
                    cellMesh( Cell).cellId = nan;
                    continue
                end
                cellMesh( Cell).cellId = double( cellList.cellId{ round}(Cell));
                
            else
                badCells = [badCells Cell];
                % Sometimes the Oufti went wrong and some cells has no information
%                     fprintf( '~~~~~~ Image %d Cell #%d has problem with its mesh! ~~~~~~\n',round, Cell)
            end
        end
        cellMesh( badCells) = [];
        cellMeshAll = [cellMeshAll; cellMesh];
    end
end

function [spotNorm, badCell] = findNormPos( pt, cellMesh, plotFlag)
    
    badCell = false;
    mesh = cellMesh.mesh;
    meshMid = cellMesh.meshMid;
    lenCum = [0; cellMesh.gridLenCum]; % add 0, so that the index matches with the mesh
    
    % calculate lNorm, match idx with mesh
    [~, dist] = findPerpFoot( pt, mesh(:, 3:4), mesh(:, 1:2)); % -: below, +: above, should be - to +
    bra = find( abs( diff( sign( dist))) == 2); % find the grid idx where distance changes sign, first & end = nan
    % dist(bra) & dist(ket) are the distance of the point to its neighboring minor axis

    % Parallel Methods
    if isempty( bra) % not sandwiched between 2 segment edges, must be in two caps then 
        if sum( dist< 0) == 0       % all positive, it's in the 1st segment
            bra = 1;
        elseif sum( dist> 0) == 0   % all negative, it's in the last segment                    
            bra = length( dist)- 1;
        else                        % the point lands on the segment!
            spotNorm = [nan nan];
            return
        end
    elseif length( bra) > 1 % bad cell, outline curls back
        badCell = true;
        spotNorm = [nan nan];
        return
    end
    ket = bra + 1;
    
    % use parallel line of segment edges to find intersection with midline
    % (LNorm) & outline (xNorm), modified @7/23/2023 by YHW 
    a = mesh( [bra ket],:);
    
    % find the mean slope of the neighboring two segments
    slopes = (a(:,4)-a(:,2))./ (a(:,3)-a(:,1));
    if sum( sign( slopes)) == 0 % two slopes have opposite sign 
        k = abs( diff( slopes))* sign( sum( slopes));
    else
        k = mean( slopes, 'omitnan');
    end
    
    [D, xPos] = findIntersect( pt, k, meshMid(ket,:), meshMid(bra,:));
    % find the intersection point of pt-D & the cell outline (cell width at lPos for this pt)
    if xPos > 0 % on the right side of the mid line
        [IntPt, ~] = findIntersect( pt, k, mesh( bra, 3:4), mesh( ket, 3:4));
    else        % on the left side of the mid line
        [IntPt, ~] = findIntersect( pt, k, mesh( bra, 1:2), mesh( ket, 1:2));
    end            
    xNorm = xPos/ norm( IntPt-D); % normalized by the width at that point                
    lPos = lenCum( bra) + norm( D- meshMid( bra,:)); % real L value by portion
    lNorm = lPos/ lenCum( end); % normalized by the total length of the major axis
    
    if abs(xNorm) > 1
        warning( 'xNorm outside the [-1 1] region, Cell %d', cellMesh.cellId)         
        xNorm = nan;
        plotFlag = true;
    end
    if abs( lNorm-0.5) > 0.5
        warning( 'lNorm outside the [0 1] region, Cell %d', cellMesh.cellId)
        lNorm = nan;
        plotFlag = true;
    end
    
%     xNorm = max( -1, min( 1, xNorm)); % set range at [-1 1]
%     lNorm = max( 0, min( 1, lNorm)); % set range at [0 1]

    spotNorm = [xNorm, lNorm];
    
    if plotFlag    
%         close
        figure()
        meshOut = cellMesh.meshOut;        
        plot( meshOut(:,1), meshOut(:,2), 'b', 'LineWidth', 1); hold on
        plot( meshMid(:,1), meshMid(:,2), 'b', 'LineWidth', 1)
        for n = 1: length( mesh)-1
            plot( mesh( n, [1 3]), mesh( n, [2 4]), 'b', 'LineWidth', 1)
        end

        % Lpos
        plot( meshMid(:,1), meshMid(:,2), 'k', 'LineWidth', 7) % total length
        plot( [ meshMid( 1:bra, 1); D(1)], [ meshMid( 1:bra, 2); D(2)],...
            'r', 'LineWidth', 3)                    
        plot( [IntPt(1) D(1)], [IntPt(2) D(2)], 'm', 'LineWidth', 3)

        % xpos
        plot( [IntPt(1) D(1)], [IntPt(2) D(2)], 'm', 'LineWidth', 10)
        plot( [D(1) pt(1)], [D(2) pt(2)], 'w', 'LineWidth', 3)                    
        scatter( IntPt(1), IntPt(2), 100, 'c', 'filled') % signal point
        scatter( D(1), D(2), 100, 'c', 'filled') % perpendicular foot
        scatter( pt(1), pt(2), 100, 'w', 'filled') % signal point
        axis equal
        figure( gcf)
        hold off
    end
%     close
end

function [D, dist] = findPerpFoot( pt, B, C)
% this function returns the coordinate of the perpendicular foot D so that 
% AD perpendicular to BC (everything in 2D) and the distance of pt to line BC
% dist > 0 if pt is on the right side of line BC (pt, B, C: counterclockwise)
% dist < 0 if pt is on the  left side of line BC (pt, B, C: clockwise)

    AB = B - pt; % vector
    BC = C - B;  % vector

    area = AB(:,1).*BC(:,2) - AB(:,2).*BC(:,1); % cross product
    side = sign( area); % -1: left side, +1: right side

    normVec = [ BC(:,2) -BC(:,1)]; % normal vector of BC in 2D 
    unitNormVec = normVec./ vecnorm( normVec, 2, 2); % unit normal vector of BC
    AD = dot(unitNormVec, AB, 2).* unitNormVec; % AD is perpendicular to BC, dot product
    D = pt + AD; % D point of intersection, Perpendicular Foot
    dist = side.* vecnorm( AD, 2, 2);
end    

function [IntPt, dist] = findIntersect( pt, k1, B, C)

    IntPt = nan( 1, 2);
    b1 = pt(2) - k1*pt(1); % pass through pt

    % kx1 + b = y1
    % kx2 + b = y2    
    % k = (y2-y1)/ (x2-x1);  b = y1 - kx1;
    k2 = ( C(2)- B(2))/ (C(1)- B(1));
    b2 = C(2) - k2* C(1);

    % k1x + b1 = y
    % k2x + b2 = y
    % x = - (b1-b2)/ (k1-k2)
    % y = k1x + b1 = k2x + b2
    IntPt(1) = -( b1-b2)/ (k1-k2);
    IntPt(2) = k2* IntPt(1) + b2;

    if sum( isnan(IntPt)) || sum( isinf(IntPt))
%         fprintf('    ~~~ Intersection has problem! ~~~\n')
        IntPt = [nan nan];
    end

    % cross product, right hand rule opposite, +1: clockwise, -1: counter        
    AB = B - pt; % vector AB
    BC = C - B;  % vector BC
    side = sign( AB(:,1).*BC(:,2) - AB(:,2).*BC(:,1)); 
    dist = side.* vecnorm( IntPt-pt, 2, 2);

end

