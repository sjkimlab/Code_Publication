clear;clc
%frameTime = 0.0217; %s

%oad( 'C:\Users\lfv\Documents\MATLAB\experimental data\SK249_tracksFinal_20220228_newAll.mat', 'tracksFinal') %change this to be the path where you have the SK249 data stored
nFrames = zeros(1, N);
for i = 1:N
    nFrames(i) = N;
end

clear tracksFinal

nFrames=nFrames(nFrames>=12);
tracksFinal = struct([]);

for ii = 1:length(nFrames)
    T_exp = nFrames(ii).*frameTime;
    params = struct('dt',dt_input,'dt_out', frameTime,'t_fin',T_exp,'l0',L+2.*R,'w0',2.*R,'totR',1,'avg',0,'D',D,'loc',1,'sigma',loc_err);
    [rne,params] = rand_diff_3D_SphCyl(params);
    tracksFinal(ii).tracksCoordXYZ = rne; 
    tracksFinal(ii).tracksCoordXY = rne(:,1:2:end);

end

tracksFinal = calc_D(tracksFinal, frameTime);


save(['D3=' num2str(D) '_D2= ' num2str(mean([tracksFinal.D12t])) '.mat'],'tracksFinal','params') %save the variables for future use

D_list = [tracksFinal.D12t];
figure;
histogram(D_list);
title('D-histogram');
xlabel('D');
ylabel('counts');
savefig('cytoplasm_D_hist.fig');
saveas(gcf,'cytoplasm_D_hist.png');

clear;clc

