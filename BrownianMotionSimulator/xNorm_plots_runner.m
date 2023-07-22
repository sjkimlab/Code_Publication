function [tracks,params] = xNorm_plots_runner(varargin)

%{
tracks contains all fields relevent to plotting x norm plots.
precent = precentage of surface particles (e.g 100 = 100% on membrane)
raw_pos = exact 3D position of particle
loc_pos = 3D position of particle after introducing localization error
xNorm = normalized x position calculated from loc_pos
h = normalized count of histogram of xNorm
mids = middle of histogram bins
ab_ratio = a/b where a is max(h) and b is h(x=0)
..._res = particles located at depth above 1/4R
..._CYL = particles located on cylindrical portion
..._res_CYL = particles located above 1/4R and is on cylindrical portion

HOW TO USE
to with defualt values:
[tracks,params] = xNorm_plots_runner()

to change any parameter: (e.g change binW to 0.01)
[tracks,params] = xNorm_plots_runner('binW',0.01)

current parameters = {'R','L','N','factor','locErr','binW','percent'}
R = (double) radius of spherical portion and width of cell
L = (double)length of cylindrical portion of cell
factor = (double) how much larger to make the mesh
locErr = (double) localization error in um
binW = (double) width of bin for histogram
percent (char) all the precentages of surface to internal particle that
will be examined.

how to plot data:
see sample plot
%}

p = inputParser;

default_R = 0.5;
default_L = 2;
default_N = 200000;
default_factor = 1.5;
default_locErr = 0.054;
default_binW = 0.05;
default_percent = ["0","10","20","30","40","50","60","70","75","80","85","90","95","100"];

addParameter(p,'R',default_R)
addParameter(p,'L',default_L)
addParameter(p,'N',default_N)
addParameter(p,'factor',default_factor)
addParameter(p,'percent',default_percent)
addParameter(p,'locErr',default_locErr)
addParameter(p,'binW',default_binW)

parse(p,varargin{:})
params = p.Results;

cc_param = struct('R',params.R,'L',params.L);
tracks = struct();

for ii = 1:length(params.percent)
    r = 0.01.*str2double(params.percent(ii));
    pos = mix(params.N,r,cc_param);
    tracks(ii).percent = params.percent(ii);
    tracks(ii).raw_pos = pos;
    XYZ = loc_error(pos,params.locErr);
    tracks(ii).loc_pos = XYZ;
    
    xnorm = xNorm(XYZ,cc_param,params.factor);
    tracks(ii).xNorm = xnorm;
    [mids,h] = myHist2(xnorm,params.binW);
    tracks(ii).h = h;
    tracks(ii).mids = mids;
    
    tracks(ii).ab_ratio = ab_ratio(h,mids);
    
end

for jj = 1:length(tracks)
    idx_list1 = tracks(jj).loc_pos(:,2)>=-0.5*cc_param.R;
    idx_list2 = abs(tracks(jj).loc_pos(:,3))<cc_param.L/2;
    idx_list3 = idx_list1 & idx_list2;
    idx_list4 = tracks(jj).loc_pos(:,2)>=-0.05*cc_param.R;
    idx_list4 = idx_list4 & idx_list2;
    
    xnorm1 = tracks(jj).xNorm(idx_list1);
    xnorm2 = tracks(jj).xNorm(idx_list2);
    xnorm3 = tracks(jj).xNorm(idx_list3);
    xnorm4 = tracks(jj).xNorm(idx_list4);
    
    [mids,tracks(jj).h_res] = myHist2(xnorm1,params.binW);
    [mids,tracks(jj).h_CYL] = myHist2(xnorm2,params.binW);
    [mids,tracks(jj).h_res_CYL] = myHist2(xnorm3,params.binW);
    [mids,tracks(jj).h_95_CYL] = myHist2(xnorm4,params.binW);
    
    tracks(jj).ab_ratio_res =  ab_ratio(tracks(jj).h_res,tracks(jj).mids);
    tracks(jj).ab_ratio_CYL =  ab_ratio(tracks(jj).h_CYL,tracks(jj).mids);
    tracks(jj).ab_ratio_res_CYL =  ab_ratio(tracks(jj).h_res_CYL,tracks(jj).mids);
    tracks(jj).ab_ratio_95_CYL =  ab_ratio(tracks(jj).h_95_CYL,tracks(jj).mids);
    
end

%{
[max,idx] = maxk(tracks(end).h(:),2);
[max1,idx1] = maxk(tracks(end).h_res(:),2);
[max2,idx2] = maxk(tracks(end).h_CYL(:),2);
[max3,idx3] = maxk(tracks(end).h_res_CYL(:),2);
tracks(14).peak_pos = (abs(tracks(end).mids(idx(1)))+abs(tracks(end).mids(idx(2))))/2;
tracks(14).peak_pos_res = (abs(tracks(end).mids(idx1(1)))+abs(tracks(end).mids(idx1(2))))/2;
tracks(14).peak_pos_CYL = (abs(tracks(end).mids(idx2(1)))+abs(tracks(end).mids(idx2(2))))/2;
tracks(14).peak_pos_res_CYL = (abs(tracks(end).mids(idx3(1)))+abs(tracks(end).mids(idx3(2))))/2;
%}

tracks(1).info = params;
    

%sample plot
figure;
%hold on
%title('surface only')
plot(tracks(end).mids,tracks(end).h);
title('xNorm')
xlabel('x-position')
ylabel('probability')
saveas(gcf,'xNorm.png')
%plot(tracks(end).mids,tracks(end).h_res);
%plot(tracks(end).mids,tracks(end).h_CYL);
%plot(tracks(end).mids,tracks(end).h_res_CYL);
%legend('whole cell','3/4 cell','CYL only','3/4 cell + CYL only')
%hold off

savefig('xNorm.fig');

end

function ratio = ab_ratio(h,mids)
a = max(h);
[minVal,idx] = min(abs(mids-0));
b = h(idx);
ratio = a/b;
end
