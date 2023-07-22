function plot_model(cc_param)
%this function will plot a graph of the spherical cylinder given cc_param 
%(a structure contaning radius info 'R', and cylinder length 'L')
%Created by Maggie Liu on 2021.8. @Kim Lab, UIUC

R = cc_param.R;
L = cc_param.L;

%define a function for the profile of capped cylinder
syms z
syms l(z)
l(z) = piecewise(-((L/2)+R)<=z<-L/2, sqrt(R.^2-(z+(L/2)).^2),...
    -L/2<=z<L/2, R,...
    L/2<=z<=((L/2)+R),sqrt(R.^2-(z-(L/2)).^2));

z = linspace(-(L/2)-R,(L/2)+R);
phi = linspace(0,2*pi);
X = double(l(z)'*cos(phi));
Z = double(l(z)'*sin(phi));
Y = repmat(z',1,100);
X = real(X);
Y = real(Y);
Z = real(Z);

%create plot
s = surf(X,Y,Z);
s.EdgeColor = 'none';
alpha(s,.3)
colormap(summer)
daspect([1 1 1])
%tracks_load = imread('tracks.png')

end