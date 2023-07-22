%The length unit is micrometer [um] and time unit is seconds [s]
clear all %;clc
%load nFramesSK187.mat %this is if we want to generate simulation tracks
%load( 'C:\Users\lfv\Documents\MATLAB\experimental data\SK249_tracksFinal_20220228_newAll.mat', 'tracksFinal') %Change this to be the location of the SK187 data file (RNase E attached to the inner membrane)

nFrames = zeros(1, N);
for i = 1:N
    nFrames(i) = N;
end

clear tracksFinal

%for our 
nFrames=nFrames(nFrames>=12);
%example 12 frames

tracksFinal = struct([])
for ii = 1:length(nFrames)
    T_exp = nFrames(ii).*frameTime;
    [track,params] = surface_random_walk('D',D,'T',T_exp,'avg',1,'R',R,'sigma', loc_err, 'L', L, 'dt_out', frameTime, 'dt', dt_input); %sigma is our localization error (um)
    tracksFinal(ii).tracksCoordXYZ = track; %with localization error
    tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
    tracksFinal(ii).pole = any(abs(track(:,3))>=L/2); %hmm... so this should be absolute value instead....
end

tracksFinal=calc_D(tracksFinal, frameTime); %this is basically calculating the D, alpha, etc from the track points, this also includes a comparison for the 3D values
%entire cell 2D
D2_withPoles = mean([tracksFinal.D12t])
D2_noPoles = mean(nonzeros(([tracksFinal.pole]-1)*-1 .* [tracksFinal.D12t]))

save(['D3=' num2str(D) '_D2= ' num2str(mean([tracksFinal.D12t])) '.mat'],'tracksFinal','params') %save the variables for future use

membrane_D = [tracksFinal.D12t]
figure;
membrane_D_hist = histogram(membrane_D);

