clear; close all;
%% PT analysis for E. col 
% %% theoretical curves
% p = [0:0.05:1]';
% a = 45; %a = L/2v using values T3' = 130s, T5' = 40s from SK98
% c = 150; %c = 1/k2 in seconds, where k2 = 0.4 per min 
% b = [c/2,c,c*2]; % b = 1/k_dPT
% est(:,1) = (a*p+b(1)*p)./(b(1)*p+c*(1-p)+a*(2-p));
% est(:,2) = (a*p+b(2)*p)./(b(2)*p+c*(1-p)+a*(2-p));
% est(:,3) = (a*p+b(3)*p)./(b(3)*p+c*(1-p)+a*(2-p));
% % diff = est(:,1)-p; diff = abs(diff); max(diff)

%% Optimizing the parameter
a = 45; %a = L/2v
c = 150; %c = 1/k2, where k2 = 0.4 per min

%% using experimental data (est)
b = c/1.9985; %This value was updated iteratively to match kpinput = kp in the end of this block
est = [27.4587260300000;56.0179174000000;47.0807768400000;0;34.4318618300000;0;1.39608109400000;0;58.9405046700000];
est = est/100; %data from steady state ratio. (1-N3/N5)
p = est*(2*a+c);
p = p./((a+b)+est*(c-b+a));

kd1 = [0.452400000000000;0.646722000000000;0.556266000000000;0.00966666700000000;0.271778000000000;0.0381396000000000;0.0304000000000000;0.0418869600000000;0.554232000000000];
%kd1 data from weak RBS strains

c = polyfit(p,kd1,1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
kd1real = c(2) %kd1* from the fit (in min-1)
kp = c(1) + c(2) %kdPT
kpinput = 1/b*60


%% error bar
b = c/1.9985;
est1 = [25.8510739900000,29.0663780800000;41.9684128900000,70.0674219100000;45.9368828100000,48.2246708800000;0,0;33.4060265400000,35.4576971200000;0,0;0,2.79216218800000;-7.70344696100000,7.70344696100000;58.2268462800000,59.6541630600000];
est1 = est1/100; %data from steady state ratio.
est = est1(:,2); %-end or +end of the error bar
p =est*(2*a+c);
p = p./((a+b)+est*(c-b+a));



% %% using experimental data (est)
% %% in case premature termination occurs right after Z5 or right before Z3
% % b = c/1.7936; %(a)
% % b = c/2.2013; %(b)
% est = [27.4587260300000;56.0179174000000;47.0807768400000;0;34.4318618300000;0;1.39608109400000;0;58.9405046700000];
% est = est/100; %data from steady state ratio.
% % case (a)
% p =est*(2*a+c);
% p = p./((b)+est*(c-b+2*a));
% 
% % % case (b)
% p =est*(2*a+c);
% p = p./((2*a+b)+est*(c-b));
% 
% kd1 = [0.452400000000000;0.646722000000000;0.556266000000000;0.00966666700000000;0.271778000000000;0.0381396000000000;0.0304000000000000;0.0418869600000000;0.554232000000000];
% %kd1 data from weak RBS strains
% 
% % c = polyfit(est,kd1,1);
% c = polyfit(p,kd1,1);
% % Display evaluated equation y = m*x + b
% disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
% kd1real = c(2)
% kp = c(1) + c(2)
% kpinput = 1/b*60