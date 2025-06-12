%{
-About-
This code contains all the work we did for the paper, Ghosal et al. (2025).
You may run each section as "selected" run.

-Author-
Aishani Ghosal, 2025 June 10
%}

%% This section is about simulation of FBM dynamics.

clear all;close all;clc;

% Setting up the parameters for the simulation
dt = 0.001; % sampling time (s)
T = 10; % total observation time (s)
frame_time  = 0.020; % camera acquisition time(s), equivalent to tE
Tstop = 0.080; % camera off time, equivalent to [gamma - tE].
%when tE = 0.02 and Tstop = 0.08, gamma = 5*tE

num_steps = T/dt;   
t = dt*(1:num_steps);  
max_lag = num_steps-1; 
lag_times = (1:max_lag) * dt;
num_intervals = floor(T / frame_time);
array = 1:num_intervals-1; 
lag_times_avg = (1:num_intervals - 1) * frame_time; % lag times for MA1 method (gamma = tE)
num_intervals2 = floor(T / (frame_time + Tstop));   
lag_times_avg2 = (Tstop + frame_time) + ((1:num_intervals2-1) - 1) *  (Tstop + frame_time); % lag times for MA2 method (gamma > tE)
X0 = log(lag_times);
X1 = log(lag_times_avg);
X2 = log(lag_times_avg2);   

% Simulation
num_trajectories = 10000;
D = 1; %um2/s^alpha
alpha = 0.25;
locError = 0.020; %um
[Y0, Y1, Y01, Y2, truetraj, trajMA, trajMALocE] = simulate_MSD(T, dt, D, alpha, locError, num_trajectories, Tstop, frame_time);    
% See simulate_MSD function for output information
 
% Simulation result plot: Figure 1b
plot_MA_only(X0, Y0, X1, Y01, t, truetraj, trajMA) ;

% Simulation result plot: Figure 1c
plot_MA_locE(X0, Y0, X1, Y01, Y1, trajMA, trajMALocE);

% Figure 2(a,b,d,e)
plot_MSD(T, X0, X1, X2, Y0, Y1, Y2, D, frame_time, alpha);


%% This section is about linear and nonlinear fitting of MSD in log-log.
% simulated data is needed.
% D and alpha are initial guesses.

% Figure 4 (a-b) - linear fitting
plot_MA_linearfitresult(Y1,Y2, X1, X2, ...
                           D, alpha);
 
% Figure 4 (c-d) - nonlinear fitting
plot_MA_nonlinearfitresult(Y1,Y2, X1, X2, ...
                          D, alpha);


%% This section is from analytical expressions for MSD' and MSD.
% Simulation is not required.

clear all;close all;clc;

% Colormaps: Figure 2(c,f) and 3(lower panel)
D = 0.01; alpha = 0.25;
optimal_param(D, alpha)
 
% colormap of the relative error of MSD in the appendix
optimal_paramspace_app(D, alpha)
 
% Figure 3 (upper panels)
frame_time = 0.02; % tE in s.
locError = 0.020; % um
plot_true_prime_MSD(D, alpha, frame_time, locError);


%% This section is about linear and noninear fitting of experimental data
% if  you want to reproduce our fig 5a

data2 = load('fig5a.dat'); % load the data file
alpha2 = 0.5;  % initial  exponent guess
D2 = 1;        % initial  diffusion coefficient guess
X2 = log(data2(:, 1)); 
Y2 = log(data2(:, 2)); 
X2 = X2';
Y2= Y2';

% nonlinear fitting parameter extraction 
[best_params21, fitted_f21] = fsimple(X2(1:floor(length(X2)/4)), Y2(1:floor(length(X2)/4)), D2, alpha2);

% linear fitting parameter extraction
[alpha_0, D_alpha_0] = estimate_diffusion_constant_and_exponent(data2(:, 2), data2(:, 1));

figure, 
De = best_params21.a;
alphae = best_params21.b;
legendText = sprintf('$$D_{e} = %.4f,\\ \\alpha_{e} = %.2f$$', De, alphae);

plot(X2, Y2,'linewidth',2, Marker='o',LineStyle='none');hold on;
plot(X2(1:floor(length(X2)/4)),fitted_f21, '--', 'linewidth',2);hold on;
ylabel('log(MSD (\mu m^{2}))','FontWeight','bold');
xlabel('log(lag time (s))','FontWeight','bold');
set(gca,'linewidth',2, 'FontSize',16,'FontWeight','bold','LineWidth',2);hold on; %xlim([-0.5 4]);hold on;
legend('1s interval',legendText,'Interpreter', 'latex','FontSize',12,'Location', 'southeast');
legend boxoff;hold on;
pbaspect([1 1 1]);


