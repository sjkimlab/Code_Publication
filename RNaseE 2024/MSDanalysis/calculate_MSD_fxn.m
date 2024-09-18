%-------------------------------------------------------------
% Author: Laura Troyer for Dr. Sangjin Kim Lab at UIUC
% Last edited date: Sept. 3, 2024
% Description: 

%======Input=======
% tracksFinal structure(from uTrack) combined from several movies
    %Note: First run format_tracksFinal
%optional: frame time (automatically 21.7 ms if no extra input is entered)
%======Output=======
% tracksFinal structure with added columns:
    %MSD: time-averaged mean squared Displacement of the individual tracks
    %Diff: apparent diffusion coefficient using the first 3 data point of the MSD
        %MSD = 4Dt + b
    %yInt: b in above equation
    %LocErr (sigma): the localization error (sigma) based on the equation 
        %MSD = 4Dt -(4/3)D*deltaT + 4*sigma^2 from https://doi.org/10.1116/1.5140087 
        %using the first 3 MSD points for fitting.
%-------------------------------------------------------------

function [tracksFinal] = calculate_MSD_fxn(tracksFinal, varargin)

    pixelSize = 160e-9; %convert localizations from pixels to SI units
    
    %get the frame time
    if isempty(varargin)
        timeStep = 21.742e-3; %units of seconds
    else
        timeStep = varargin{1}; %units of seconds
    end
    
    %calculate MSD, D, and error for each track in tracksFinal
    nTracks = length(tracksFinal);
    for i = 1:nTracks
        %get the track length for this track
        if isfield(tracksFinal, 'tracksCoordAmpCG')
            nFrames = size(tracksFinal(i).tracksCoordAmpCG,2)/8;
        else
            nFrames = size(tracksFinal(i).tracksCoordXY ,1); %if tracksCoordAmpCG was not kept in the tracksFinal data set
        end

        %load in the track position data 
        if isfield(tracksFinal, 'tracksCoordAmpCG')
            Traj = zeros(nFrames,2); % to store x and y coordinates
            Coords=tracksFinal(i).tracksCoordAmpCG;
            Traj(:, 1) = Coords(1:8:end) * pixelSize;
            Traj(:, 2) = Coords(2:8:end) * pixelSize;
        else
            Traj = tracksFinal(i).tracksCoordXY;
        end
        
        %calculate this track's MSD
        nTau = nFrames-1; %max lag time for this track
        msd = zeros(nTau,1);

        for tau = 1:nTau
            dr = Traj(1+tau:end, 1:2) - Traj(1:end-tau, 1:2);
            dr2 = sum(dr.^2, 2);
            msd(tau)=mean(dr2); %units of m^2 
        end
        
        %input the MSD value into the tracksFinal structure
        [tracksFinal(i).MSD] = msd;

        %calculate this track's D coefficient. Since we calculate D
        %from the first 3 MSD points and want to using only points in
        %the first 25% of the track: only calculate D if the track is
        %at least 12 frames long
            
        if nFrames >= 12 %check the track is at least 12 frames long                
            %calculate D and localization error for this track
            [fitobject,gof,output] = fit(((1:3)*timeStep)', (msd(1:3)), 'poly1');
            MyCoeffs = coeffvalues(fitobject); %MyCoeffs is (1) the slope and (2) the y-intercept from linear fit
            D = MyCoeffs(1)/4; %D value from 4Dt term, units of m^2/s
            yInt = MyCoeffs(2); %y-Intercept from linear fit. Includes dynamic and static error, units of m^2
            LocErr= sqrt(abs(MyCoeffs(2)/4)); %(static) localization error sigma, units of m

            %insert these values into the tracksFinal structure
            tracksFinal(i).D = D;
            tracksFinal(i).yInt = yInt;
            tracksFinal(i).LocErr = LocErr;
            
        end
    end
end