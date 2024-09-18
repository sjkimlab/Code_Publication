function [newX, Comb] = num_xNorm_fast( wid, fCut, locErrX, locErrY, dilF, MBper, approxFlag)
%{
Author: Yu-Huan Wang   (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 4/13/2023
    Last update date: 10/10/2023

Description: this script model the circular cross-section of a cell as a
circle of radius 1, x & y ranges from -1 to 1, all operations (locErr
fCut) are performed on this circle but with a scaling factor that takes the
cell width & dilF into consideration, the xNorm was eventually scaled again
to match the real data.

Note: In the paper, the circular cross-section is in the x-z plane. Here, I
used y instead of z

------ Input ------
wid:        Cell width (um)
fCut:       Defocused region at the bottom part of the cell (um)
locErrX:    Localization error in x axis
locErrY:    Localization error in y axis
dilF:       Dilation factor, relative difference between inner membrane vs. cell outline
MBper:      Membrane binding percentage (MB%)
approxFlag: Flag for approximating the integral function

------ Output ------
newX:       Bin centers for numerical xNorm histogram
Comb:       Corresponding xNorm histogram value for each bin


10.10: added approxFlag feature, doesn't perform integral every time
---------------------------------------------------------------------------
%}    
    
    global funY0
    
    % example of parameters
%     wid = 1; dilF = 1.6; MBper = 50; locErrX = 100; locErrY = 100; fCut = 0.3;
    
    % binWidth is 2/nBin, x ranges from [-1,1]
    nBin = 100;    binEdges = linspace( -1, 1, nBin+1); % assumed r=1
    tmp = movmean( binEdges, 2);    binCenters = tmp(2:end); 
    
    % [Membrane]     calculate the surface density (arc length) in each bin, normalized later
    surf = -diff( acos( binEdges)); % arc = r* d(theta), theta = arccos( x/r), r=1
    
    % [Cytoplasmic]  calculate the area in each bin, normalized later
    y = sqrt( 1- binCenters.^2); % corresponding y of each bin centers
    cyto = y;
    
    % rescale the model cell to match real inner membrane,  r=1 --> R/dilF
    R = wid/2;  scaleF = dilF/ R; % scaleF ~ 3
    
    % fCut is determined from the bottom of the cell ~ 0.3 um
    % fCutHere is relative to the cell midline, with 0 at the middle of cell
    fCutHere = (R- fCut)* scaleF; % for integral later
    
    % Gaussian blue vector [1,51] in x axis, the locErr was scaled by scaleF
    blurX = normpdf( binEdges, 0, locErrX* scaleF/ 1000);
    
    % cdf of the blur function in y axis
    funY = @(x) normcdf( x, 0, locErrY* scaleF/ 1000);
    
    % [Membrane Part Blur]
    surfCut = surf.* ( funY( fCutHere + y) + funY( fCutHere - y))/2; % bottom + top part
    surfFinal = conv( surfCut, blurX); % apply locErr in x axis
    
    % [Cytoplasmic Part Blur]
    if approxFlag % interpolate the integral function to speed up calculation
        if isempty( funY0)
            funY0 = getApproxFunY(); % define interpolation function
        end
        sigma = locErrY* scaleF/ 1000;
        cytoPart = (funY0( (fCutHere+y)/ sigma) - funY0( (fCutHere-y)/ sigma))* sigma./ (2*y); % change of variables
    else
        % performe integral for each bin
        cytoPart = nan(1, length( y));    
        for i = 1: ceil( nBin/2)
            cytoPart( i) = integral( funY, fCutHere-y(i), fCutHere+y(i))/ (2*y(i));
        end
        cytoPart( i+1: end) = flip( cytoPart( 1: floor( nBin/2))); % use symmetry to speed up
    end
    
    cytoCut = cyto.* cytoPart; % apply fCut with locErrY
    cytoFinal = conv( cytoCut, blurX); % apply locErr in x axis
    
    % conv: [1,50]*[1,51]-->[1,100], new binEdges:[-2,2] with 2*nBin
    conX = linspace( binCenters(1)-1, binCenters(end)+1, nBin*2); % x coordinate after convolution
    
    % Combine two parts based on the ratio before applying fCut [surf,cyto], normalization
    %     sum( conv(a,blur)) = sum(a)* sum(blur)
    %     surf/sum(surf)  vs.  cyto/sum(cyto)  two parts contributes equal
    Comb = ( surfFinal/ sum( surf)* MBper + cytoFinal/ sum( cyto)*( 100-MBper))/ 100;
    
    newX = conX/ dilF; % real bin centers after scaling by dilF (membrane/outline)
    
end

function funY0 = getApproxFunY()
    r = [-50:-4 -3:0.001:3 4:50];
    
    y = nan( 1, length( r));
    for i = 1: length( r)
        y( i) = integral( @(x) normcdf(x, 0, 1), 0, r(i));
        % normcdf(x,mu,sigma) = normcdf((x–mu)/sigma,0,1), they are equivalent
    end
    
    funY0 = @(x) interp1( r, y, x);
end
