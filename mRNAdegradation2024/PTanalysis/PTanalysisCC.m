%% PT analysis for C. crescentus (SK438)
close all; clear;

m = [0:0.1:10]/60; % candidate for kd1 (in sec-1)
est = 0.16; % equals to 1-N3/N5, estimation of PT based on steady-state level difference
t = 75; %L/2v = 150s/2
kd2 = 0.4/60; % of 3' [unit: sec-1]


% kd1 = 1.28/60
% kdPT = 1.28/60; % of 5'
% a = exp(-kd1*t);
% PT = (est-1)*(t*a + t*a*a)+1/kd2;
% PTr = (est-1)*(t*a*a-a/kdPT)+1/kd2;
% PT = PT/PTr;
% 
% %1.28/60 = kd1*(1-PT)+kdPT*PT;
% PT
% kd1out = (1.28/60-kdPT*PT)/(1-PT)


%---parameter scan
for i = 1:length(m)
    kd1= m(i);
    for j = 1:length(m)
        kdPT = m(j);
        
        a = exp(-kd1*t);
        PT = (est-1)*(t*a + t*a*a)+1/kd2;
        PTr = (est-1)*(t*a*a-a/kdPT)+1/kd2;
        PTm(i,j) = PT/PTr; %expected PT
        kd1m(i,j) = (1.28/60-kdPT*PT)/(1-PT); %another kd1 estimation based on equation (1)
        kd1e(i,j) = (1.28/60-kdPT*PT)/(1-PT)-kd1; %difference
    end;
end;
        
M = abs(kd1e); [R,C] = find(M==min(M(:)))
kd1e(R,C), kd1m(R,C), PTm(R,C)
