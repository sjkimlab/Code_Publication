function [output_track,param] = surface_random_walk(varargin)
% This function will simulate the random walk process on a cell modeled by a
% cylinder with spherical caps with user specified diffusion parameters.
% The magnitude of each step is drawn from rayleigh distribution.
%
% INPUT:
% varargin - this is a variable-length input argument list. This is used to
% specify settings related to the simulation. This include geometrical description
% of the cell, simulation time steps and diffusion coefficient.
% See code body for a full description of available options.
%
% OUTPUT :
% output_track - this is the path taken by the walker. The first row is the
% XYZ position of the walker at the first output time step, the second row
% is the XYZ position at the second output, and so on.
% param - this is the input parameters specified by the user.
%
% Units: length measurement are in micrometers [um], time measurement are
% in seconds [s]
%
% Created by Maggie Liu on 2021.9. @Kim Lab, UIUC
% Last revised 2022.3.24
%% parse input
p = inputParser;

%radius of sphere
default_R = 0.5;

%length of cylinder *total cell length is L+2R*
default_L = 2; 

%"silent" steps used for finner grained simulation. positions at the silent
%steps will not become an output. 
default_dt = 0.007; %motion bluring timing

%output time step. %This is the frame rate of the experiment.
default_dt_out = 0.021;

%amount of time that the camera is exposed to light.  This should be a multiple of 'dt' and should match the
%frame speed of experimental aperatus.
default_exposureTime = 0.0217;

%total simulation time. 
default_T = 0.4;

%static localization error
default_sigma = 0;

%input diffusion coefficient
default_D = 0.2;

%if 'avg' is true, then positions visited by the walker in the silent cycle
%will be averaged and recorded as output position. if 'avg' is false, then
%only the last position visited by the walker in the silent cycle will be 
%recorded as output position
default_avg = true;

%if 'origin' is true, then every track will begin at (0,R,0) (this
%traslates to projected xy coordinate of (0,0). Setting this field as true
%with a very large radius (e.g R=100) can simulate the case where particles
%are diffusing on a infinite plane) 
default_origin = false;



positive = @(x) isnumeric(x) && (x >= 0);
isBoo = @(x) x==1 || x==0;

addParameter(p,'dt',default_dt,positive)
addParameter(p,'dt_out',default_dt_out,positive)
addParameter(p,'T',default_T,positive)
addParameter(p,'sigma',default_sigma,positive)
addParameter(p,'D',default_D,positive)
addParameter(p,'avg',default_avg,isBoo)
addParameter(p,'R',default_R)
addParameter(p,'L',default_L)
addParameter(p,'origin',default_origin)
addParameter(p,'exposureTime', default_exposureTime)

parse(p,varargin{:})

param = p.Results;

%geometriic parameters
R = param.R; %radius of cell
L = param.L; %length of cylinderical part (total cell length is L+2R)

%load parameters 
D = param.D;
dt_out = param.dt_out;
dt = param.dt;
T = param.T;
sigma = param.sigma;
avg = param.avg;
origin = param.origin;
exposureTime = param.exposureTime;

%%define some convineint variables 
tot_steps = floor(T/dt);
output_steps = floor(T/dt_out)+1;
tot_track = zeros(tot_steps,3);
output_track = zeros(tot_steps,3);
cc_param = struct('R',R,'L',L);

%define my pdf function from which step sizes are drawn
px = linspace(0,0.5);
p = px./(2.*D.*dt).*exp(-px.^2./(4*D*dt));
if D > 0
    step_size=randpdf(p,px,[tot_steps,1]);
else
    step_size=zeros(tot_steps,1);
end

cycle = floor(dt_out/dt); %number of silent cycles per output
pos_list = zeros(cycle,3); %same as above, but for z values

%pick random starting point
if origin
    pos = [0,R,0];
else
    pos = surface_pt(cc_param,1);
end
tot_track(1,:) = pos;
output_track(1,:)=pos;
pos_list(1,:) = pos;

rng('shuffle');

for ii = 2:tot_steps
    
    %get initial point position info
    z0 = pos(3);
    phi0 = atan2(pos(2),pos(1));
    
    rand_angle = 2*pi*rand;
    
    %if initial point is on the sphere
    if abs(z0)>(L/2)
        if z0>(L/2) %upper sphere
            theta0 = atan2(l(z0),(z0-(L/2)));
        else %lower sphere
            theta0 = atan2(l(z0),(z0+(L/2)));
        end
        
        %see Rodrigues' rotation formula. the point of below is just to
        %locate a point with set distance to initial point
        phi_perp = phi0+pi./2;
        rand_pos_sph = [R,rand_angle,step_size(ii)/R];
        [x,y,z] = sph2rec(rand_pos_sph);
        v = [x,y,z];
        k = [cos(phi_perp),sin(phi_perp),0];
        v_rot = v.*cos(theta0) +cross(k,v).*sin(theta0) + k.*(sum(k.*v)).*(1-cos(theta0));
        
        if z0>(L/2) %add back the shift
            pos = v_rot + [0,0,L/2];
        else
            pos = v_rot + [0,0,-L/2];
        end
        
        %if the walk ended up on the cylinder, I correct for the radius
        if abs(pos(3))<= L/2
            phi_temp = atan2(pos(2),pos(1));
            pos = [R.*cos(phi_temp),R*sin(phi_temp),pos(3)];
        end
        
    end
    
    %if initial point is on the cylinder
    if abs(z0)<=(L/2)
        step_z = step_size(ii)*cos(rand_angle);
        step_phi = step_size(ii)*sin(rand_angle);
        pos = [R.*cos(phi0+step_phi./R),R.*sin(phi0+step_phi./R),z0+step_z];
        if abs(step_z+z0)>(L/2) %if point end up on sphere I correct for radius
            pos = [l(z0+step_z).*cos(phi0+step_phi./R),l(z0+step_z)*sin(phi0+step_phi./R),z0+step_z];
        end
    end
    
    %update track
    tot_track(ii,:) = pos;
    
    %code below codes for position averaging aka motion blurring
    idx = mod(ii-1,cycle) + 1;
    pos_list(idx,:) = pos;
    lastPoint = round(exposureTime/dt);
    if mod(ii,cycle)== 1
        if avg
            
            r = (mean(pos_list(1:lastPoint,3)));
            p = l(r)*(sum(pos_list(1:lastPoint,1)))./sqrt((sum(pos_list(1:lastPoint,1)))^2+(sum(pos_list(1:lastPoint,2)))^2);
            q = (sum(pos_list(1:lastPoint,2)))*p/(sum(pos_list(1:lastPoint,1)));

            output_track(ii,:) = [p,q,mean(pos_list(1:lastPoint,3))];
            
        else
            output_track(ii,:) = pos;
        end
    end
end

%remove zero elements of output track
output_track( ~any(output_track,2), : ) = [];
%output_track = output_track(1:output_steps,:);

%localization error
output_track = output_track + sigma.*normrnd(0,1,[size(output_track,1),3]);

    function [x,y,z] = sph2rec(pos_sph) %converts spherical coord to cartesian coord
        phi = pos_sph(2);
        theta = pos_sph(3);
        x = R.*cos(phi).*sin(theta);
        y = R.*sin(phi).*sin(theta);
        z = R.*cos(theta);
    end

    function r = l(z) %return radius of capped cylinder
        if z>=L/2 && z<=(L/2+R)
            r = sqrt(R.^2-(z-L/2).^2);
        elseif -L/2>=z && -(L/2+R)<=z
            r = sqrt(R.^2-(z+L/2).^2);
        elseif -L/2<=z && z<=L/2
            r = R;
        else
            r= 0;
        end
    end
%close;

%if nargout==0 %if no output requested a 3D plot will be generated for visualization
%    figure(1)
%    plot_model(cc_param);
%    hold on
%    plot3(tot_track(:,1),tot_track(:,2),tot_track(:,3),'-')
%    daspect([1 1 1])
%    plot3(output_track(:,1),output_track(:,2),output_track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)
%    daspect([1 1 1])
%    set(gcf,'color','w');
%    set(gca,'linewidth',2)
%    savefig('tracks.fig');
%    saveas(gcf,'tracks.png');
    %tracks_load = imread('tracks.png');
%end

end

function rand_pt = surface_pt(cc_param,N) %generates a random point on the surface

R = cc_param.R;
L = cc_param.L;
shift = (L/2);

rand_num = unifrnd(-(shift+R),(shift+R),N,1);
rand_angle = unifrnd(0,2*pi,N,1);
rand_pt = zeros(N,3);

for ii = 1:N
    r = l(rand_num(ii));
    rand_pt(ii,:) = [r*cos(rand_angle(ii)),r*sin(rand_angle(ii)),rand_num(ii)];
end

    function r = l(z)
        if abs(z)<= L/2
            r = R;
            return
        end
        
        if z>L/2
            r = sqrt(R.^2-(z-L/2).^2);
            return
        end
        
        if z<-L/2
            r = sqrt(R.^2-(z+L/2).^2);
            return
        end    
    end

end

function x=randpdf(p,px,dim)
% Copyright (c) 2009, Adam Nieslony
% All rights reserved.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% RANDPDF
%   Random numbers from a user defined distribution
%
% SYNTAX:
%   x = randpdf(p, px, dim)
%       randpdf(p, px, dim)
%
% INPUT:
%   p   - probability density,
%   px  - values for probability density,
%   dim - dimension for the output matrix.
%
% OUTPUT:
%   x   - random numbers. Run function without output for some plots.
%
% DESCRIPTION:
%   x = randpdf(p, px, dim) returns the matrix of random numbers from
%   probability density distribution defined in p and px. p are the density
%   (the y axis) and px are the value (the x axis) of the pdf. p and px
%   must be of the same length.
%   dim define the output matrix dimensions, for example dim=[100 3] define
%   the 100x3 two dimensional matrix with 300 random numbers.
%
%   REMEMBER: This is not a realy random number generator but only
%   some kind of transformation of uniformly distributed pseudorandom
%   numbers to desired pdf!
%
% By Adam Nies³ony, Opole University of Technology, Poland

% check the number of input
error(nargchk(3, 3, nargin))

% vectorization and normalization of the input pdf
px=px(:);
p=p(:)./trapz(px,p(:));

% interpolation of the input pdf for better integration
% in my opinion 10000 point is sufficient...
pxi=[linspace(min(px),max(px),10000)]';
pi=interp1(px,p,pxi,'linear');

% computing the cumulative distribution function for input pdf
cdfp = cumtrapz(pxi,pi);

% finding the parts of cdf parallel to the X axis
ind=[true; not(diff(cdfp)==0)];

% and cut out the parts
cdfp=cdfp(ind);
pi=pi(ind);
pxi=pxi(ind);

% generating the uniform distributed random numbers
uniformDistNum=rand(dim);

% and distributing the numbers using cdf from input pdf
userDistNum=interp1(cdfp,pxi,uniformDistNum(:)','linear');

x=reshape(userDistNum,dim);
end

%function plot_model(cc_param)
%this function will plot a graph of the spherical cylinder given cc_param 
%(a structure contaning radius info 'R', and cylinder length 'L')
%Created by Maggie Liu on 2021.8. @Kim Lab, UIUC

%R = cc_param.R;
%L = cc_param.L;

%define a function for the profile of capped cylinder
%syms z
%syms l(z)
%l(z) = piecewise(-((L/2)+R)<=z<-L/2, sqrt(R.^2-(z+(L/2)).^2),...
%    -L/2<=z<L/2, R,...
%    L/2<=z<=((L/2)+R),sqrt(R.^2-(z-(L/2)).^2));

%z = linspace(-(L/2)-R,(L/2)+R);
%phi = linspace(0,2*pi);
%X = double(l(z)'*cos(phi));
%Y = double(l(z)'*sin(phi));
%Z = repmat(z',1,100);

%create plot
%s = surf(X,Y,Z);
%s.EdgeColor = 'none';
%alpha(s,.3)
%colormap(summer)
%daspect([1 1 1])

%tracks_load = imread('tracks.png')

%end

