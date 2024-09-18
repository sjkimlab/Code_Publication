function [eaTaMSD, tau, eaTaMSD_SD, eaTaMSD_SEM, taMSD] = calc_eaTaMSD(tracksFinal, timeStep)
    %This is a function to create a matrix of TA MSD from multiple tracks
    %in tracksFinal and then perform an ensemble average to get the EATA
    %MSD. This script was last updated by Laura Troyer for use in Dr.
    %Sangjin Kim's lab at UIUC on Sept. 5, 2024.
    
    %==========================================
    %Input: 
    %tracksFinal structure. (Output from Oufti and after calculating MSD
        %values. LT uses calculate_MSD_fxn.m. 
    %timeStep: This is the frame time (aka lag time, tau) for the data set,
    %should be a 1x1 double
    %============================================
    %Output:
    %eaTaMSD: 1x50 vector of the eaTA MSD data. The ith item corresponds to
        %the eaTAMSD value at the i*timeStep time apart.
    %tau: 1x50 vector of the lag times at the i*timeStep time apart.
    %eaTaMSD_SD: 1x50 vector of the standard deviation of the TA MSD values
        %at the i*timeStep time apart
    %eaTaMSD_SEM: 1x50 vector of the standard error of the mean of the TA
        %MSD values at the i*timeStep time apart
    %taMSd: nx50 matrix where n is the number of tracks in the tracksFinal
        %structure that are at least 12 frames long. Note that for any track
        %that has less than 51 frames, NaN values will be put as placeholders
        %where there are no taMSD values for that specific track. For any
    %tracks longer than 51 frames, only the first 50 taMSD values will be
        %recorded.
    %==============================================
    
    taMSD = []; % to store TA MSD
   
    for i = 1:length(tracksFinal)
        %calculate the length of the ith track
        if isfield(tracksFinal, 'tracksCoordXY')
            trackLen = length(tracksFinal(i).tracksCoordXY);
        else
            trackLen = length(tracksFinal(i).tracksCoordAmpCG)/8;
        end


        if trackLen >= 12 %set minimum track length to include in taMSD.
            new = nan(50,1); %create a vector for this track's taMSD
            if trackLen < 51 %these frames have less than 50 lag times for MSD, so will need to include NaN values as placeholders for larger lag times
                new(1:trackLen-1) = tracksFinal(i).MSD*1e12;
            else %these tracks have MSD values for the entire 50 lag times  
                new = tracksFinal(i).MSD*1e12;
                new = new(1:50);
            end
            [r, c] = size(new);
            if r>1 %this if statement accounts for if the tracksFinal(i).MSD values are in a column vecotor or row vector and ensures each new vector is oriented in the same direction
                taMSD = [taMSD; new'];
            else
                taMSd = [taMSD; new];
            end
        end
    end

    %calcualte the EATA MSD values along with the standard deviation and
    %the standard error of the mean for each EATA MSD lag time.
    eaTaMSD = nanmean(taMSD);
    tau = [1:50] * timeStep;
    eaTaMSD_SD = nanstd(taMSD);
    eaTaMSD_SEM = nanstd(taMSD)./sqrt(sum(~isnan(taMSD)));
end