%to calculate D, alpha, and MSD of tracks
function tracksFinal = calc_D(tracksFinal, frameTime, twenty_five_percent)

for ii = 1:length(tracksFinal)
track = tracksFinal(ii).tracksCoordXY;
[MSD,~] = MSD_2D(track);
tracksFinal(ii).MSD = MSD;
t_axis = frameTime:frameTime:length(track)*frameTime;

[tracksFinal(ii).D,tracksFinal(ii).LocEr,tracksFinal(ii).Alpha] = Diff(t_axis,tracksFinal(ii).MSD,twenty_five_percent);
end

function [D,locErr,alpha] = Diff(t_axis,MSD,t)
    [slope] = polyfit(t_axis(1:t),MSD(1:t),1);
    D = slope(1)/4;
    locErr = sqrt(abs(slope(2)/4));
    [slope] = polyfit(log(t_axis(1:t)),log(MSD(1:t)),1);
    alpha = slope(1);
end

function [AVG,MSD_matrix] = MSD_2D(track)
len = length(track);
MSD_matrix = zeros(len-1);
for i = 1:len-1
    A = track(i,:);
    for jj = 1:len-i
        B = track(jj+i,:);
        MSD_matrix(i,jj) = sum((A-B).^2);
    end
end
AVG = zeros(len-1,1);
for i = 1:len-1
AVG(i) = mean(nonzeros(MSD_matrix(:,i)));
end
end
end