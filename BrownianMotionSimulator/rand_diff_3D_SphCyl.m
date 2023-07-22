function [rne,params] = rand_diff_3D_SphCyl(params)
% To simulate 3D dynamics of RNaseE in 249 and 249-rif
% This is a simplified version of the code used in the paper
% Surovtsev et al. PNAS 2016 113 E7266 
% http://www.pnas.org/content/113/46/E7268.full
%
% Simulations organized as a nested loops with "silent" internal loop that run simulations without saving coordinates, 
% and external loop in which corrdinates are saved at each iteration.
% All simulations are within a 3D spherocylinder with reflective or sticky 
% boundaries (default is reflective).
%
% This simulation is meant to mimic the measurement
% process. Therefore the reported molecule positions after every frame are
% taken as the average of all the microstep positions, and dynamic 
% localization error is modeled as added gaussian noise. In this sense the 
% reported positions are like the "measured" positions from under a
% microscope. (params.avg = true, params.loc = true)
%
%
% INPUT:
%   params - parameters for simulation with multiple fields, for example:
%     params.dt - Simulation time step
%     params.dt_out - data output time step
%     params.t_fin - end time of simulation
%   Full list of parameters can be found in a nested function "change_parameters" below
%   Default values will be used for parameters value/fields not provided in input
%
% OUTPUT:
%   rne - measured coordinates of RNaseE vs time, i-th row is a snap shot of coordinates of all RNaseE (X1,Y1, Z1, X2,Y2, Z2 ...) at a given time,(i-1)*dt_out 
%   params - parameters used for simulation, also contain version/date info
%   NOTE: When getting/using data from simulation, discard position at time 
%   t = 0, because it takes into account neither averaging due to finite camera
%   exposure nor localization error.
%   
% Andrew Maytin
% 8.2.2020
%__________________________________________________________
%%


%% Parameters
 % First, define all required parameters of simulations 
 % If values provided in the input "params" - change default values to ones provided in params
 % Note: parameters units should be self-consistent. For example, use seconds and micrometers as all time and length units  
 script=mfilename; % save code name for the output 
 params = change_parameters(params); % get all parameters values

% Split 'params' into individual parameters
 % simulation time parameters
dt=params.dt;  % Simulation time step
dt_out=params.dt_out; % data output time step
t_fin=params.t_fin; % end time of simulation
 % simulation cell parameters
l0=params.l0;    % cell length
w0=params.w0;  % cell width

totR=params.totR;  % number of RNaseE
D=params.D; % Diffusion coefficient
loc=params.loc;
avg=params.avg;
sigma=params.sigma; % dynamic localization error
boundary= params.boundary;
reflective = params.reflective;

%% INTIALIZATION
rng('shuffle'); % to randomly reset RND-generator and change generator method from default

% geometric parameters
l00=l0/2; % long-axis
w00=w0/2;% short axis 
Rcyl=w00; % Radius of the cyliner and caps
Lcyl0=l00-Rcyl; % half-length of cylinder part
 
% setting initial coordinates for RNaseE, randomly distributed   
rand_XYZ=get_rand_XYZ(totR); % get random X, Y and Z-s
x_R0=rand_XYZ(1,:);  
y_R0=rand_XYZ(2,:); 
z_R0=rand_XYZ(3,:); 

x_R=x_R0; 
y_R=y_R0; 
z_R=z_R0; 
 
% making data collectors for the output
rne=111*ones(floor(t_fin/dt_out)+1,3*totR);  % RNaseE coordinates

% saving initial values - RNaseE coordinates states - at t=0;  
% NOTE: When getting/using data from simulation, discard position at time 
% t = 0, because it takes into account neither averaging due to finite camera
% exposure nor localization error.
rne(1,1:3:end)=x_R;
rne(1,2:3:end)=y_R;
rne(1,3:3:end)=z_R;

%% SIMULATION

ii0=1; % counter for data saving

%pause on

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
for tt1=dt_out:dt_out:t_fin % "save cycle", save data at each step
 
  % Get pool of random steps, i.e., Brownian dynamics steps
  % Used the formula from https://www.nature.com/articles/nmeth.2367
  Dpool=sqrt(2*D*dt)*normrnd(0,1,[1,3*totR*ceil(dt_out/dt)]);  

  % counter for drawing numbers from the pool
  ddB=0; 
  
  % arrays to store the microtrajectories
  x_M = ones(floor(dt_out/dt),totR);
  y_M = ones(floor(dt_out/dt),totR);
  z_M = ones(floor(dt_out)/dt,totR);
  
  % counter for storing microtrajectories
  ii1 = 0;
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for tt2=dt:dt:dt_out % "silent" cycle, move without saving data
      
    % ~~~~~ Make BD moves    
 
    % - RNaseE moves
       xB_new=x_R+Dpool(1,ddB+1:ddB+totR);
       yB_new=y_R+Dpool(1,ddB+totR+1:ddB+2*totR);
       zB_new=z_R+Dpool(1,ddB+2*totR+1:ddB+3*totR);
       if boundary
       % apply reflecting or sticky boundaries (default: reflective)
       [x_R,y_R,z_R] = apply_boundaries(xB_new,yB_new,zB_new,reflective);
       else
       x_R=xB_new;
       y_R=yB_new;
       z_R=zB_new;
       end
       
       ii1 = ii1+1;
       x_M(ii1,:) = x_R;
       y_M(ii1,:) = y_R;
       z_M(ii1,:) = z_R;
     
     % shift the pool counter 
     ddB=ddB+3*totR;
           
  end % "silent" cycle
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   % Current Data collection
   ii0=ii0+1;
   
   if avg
   % position saved as the centroid of microtrajectories to mimic analysis procedure of experimental images
   rne(ii0,1:3:end)=mean(x_M);
   rne(ii0,2:3:end)=mean(y_M);
   rne(ii0,3:3:end)=mean(z_M);
   else
   % position saved as final position after all microsteps
   rne(ii0,1:3:end)=x_R;
   rne(ii0,2:3:end)=y_R;
   rne(ii0,3:3:end)=z_R;
   end
   
   if loc
   % dynamic localization error is applied to each centroid location in both x and y coordinates 
   rne(ii0,1:3:end)=rne(ii0,1:3:end)+sigma*normrnd(0,1,[1,totR]);
   rne(ii0,2:3:end)=rne(ii0,2:3:end)+sigma*normrnd(0,1,[1,totR]);
   rne(ii0,3:3:end)=rne(ii0,3:3:end);
   end

end % "save" cycle
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%% NESTED FUNCTIONS
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function rand_XYZ=get_rand_XYZ(nn)
% returns n random XYZ triplets ([x1,x2,x3...; y1,y2,y3...; z1,z2,z3...]) within spherocylinder
 uni_rand=rand(3,(3+(nn<100)*100)*nn);
 rand_0x=l0*uni_rand(1,:)-l00;
 rand_0y=w0*uni_rand(2,:)-w00;
 rand_0z=w0*uni_rand(3,:)-w00;
 xL=abs(rand_0x)-Lcyl0;
 r_old=sqrt(heaviside(xL).*xL.^2+rand_0y.^2+rand_0z.^2);
 rand_1x=rand_0x(r_old<Rcyl);
 rand_1y=rand_0y(r_old<Rcyl);
 rand_1z=rand_0z(r_old<Rcyl);
 rand_XYZ(1,:)=rand_1x(1:nn);
 rand_XYZ(2,:)=rand_1y(1:nn);
 rand_XYZ(3,:)=rand_1z(1:nn);

end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function params_out = change_parameters(params_in)
% declare default values of parameters and change them to values specified in params_in 
  % declare default parameters values
 params0.dt=0.1;  % Simulation time step
 params0.dt_out=1; % data output time step
 params0.t_fin=10; % end time of simulation

 l_0=3.0; w_0=1.1; % cell tip-to-tip length and cell width
 params0.l0=l_0; %l_0=2.5;    % cell length
 params0.w0=w_0; %w_0=0.5;  % cell width

 params0.totR=90;  % number of RNaseE

 params0.D=0.2; % diffusion coefficient
 
 params0.boundary=false; %boolean to activate boundary
 params0.reflective=false; %boolean to activate boundary reflection
 
 params0.avg=false; %boolean to activate averaging
 
 params0.loc=false; %boolean to activate dynamic localization error
 params0.sigma=0.00; % dynamic localization error
  
 params0.space='um'; % space units
 params0.time='s'; % time units
 params0.dim=3; % Dimension of simulations;
 params0.script=script; %'random_diffusion_3D_SphCyl';
 params0.model='1-state model';
 params0.version='2020.08.02';
  
 % check what is provided as input parameters, i.e., in 'params_in', and substitute those in default values
 fldnms_in = fieldnames(params_in);
 fldnms_0 = fieldnames(params0);
 params_out=params0; 
 for kk = 1:length(fldnms_in)
    not_found=1;
    for jj = 1:length(fldnms_0)
      if (strcmpi(fldnms_in{kk},fldnms_0{jj}))
         params_out.(fldnms_0{jj}) = params_in.(fldnms_in{kk});
         not_found=0;
      end
    end
    if not_found
        disp(['Warning: parameter ',fldnms_in{kk},' not found in the list of parameters, check spelling...'])
    end
 end
 if params_out.dt_out<params_out.dt
    params_out.dt_out=params_out.dt; 
 end
  
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function [x_new,y_new,z_new] = apply_boundaries(x,y,z,reflective)  
% applies reflective (or sticky) boundaries shaped as spherocylinder with cylinder length Lcyl and cylinder and caps radious Rcyl
   x_new=x;
   y_new=y;
   z_new=z;

  xL=abs(x)-Lcyl0;
  r_old=sqrt((xL>0).*xL.^2+y.^2+z.^2); % if x is in spherical cap, calculate distance from center of sphere
                                       % if x is in cylinder, calculate distance from centerline of cylinder                                                                      
  ind=(r_old>Rcyl);  
  if reflective
  r_new=2*Rcyl-r_old(ind);
  else                                 % if reflective is false, boundary is "sticky": particle position becomes 
  r_new=Rcyl;                          % stuck to boundary for that microstep.
  end  
  x_new(ind)=(xL(ind)>0).*sign(x(ind)).*(Lcyl0+xL(ind).*r_new./r_old(ind)) + (xL(ind)<=0).*x(ind);
  y_new(ind)=y(ind).*r_new./r_old(ind);
  z_new(ind)=z(ind).*r_new./r_old(ind);
end

%% END of ENDS

end