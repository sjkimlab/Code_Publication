function Xans = getAffineFactor(YFPfilepath, Cy5filepath, Cy3filepath)
%{
-About-
This code Affine transformation factor to correct for (X,Y) shift between Cy5 and Cy3 channel.
Tetraspek beads were imaged in YFP, Cy5, and Cy3 channels.
Then, the images were analyzed by spotFinderF of microbeTracker.
spotFinderF identifies spots, in the absence of cell mesh information.
spotList were saved as YFPspot.mat, Cy5spot.mat, Cy3spot.mat.

One may need to change "pixel" variable, "upperlimit" variable and "threshold" variale.

-Inputs-
YFPfilepath, Cy5filepath, Cy3filepath

-varargin-

-Outputs-
Xans: Affine transformation factor

-Example-
Xans = getAffineFactor('C.Beadcontrol\YFPspot','C.Beadcontrol\Cy5spot','C.Beadcontrol\Cy3spot')
see FISHworkflow for the input
  
-Supplementary-

-Keywords-
two-color, colocalization, channel correaction

-Dependencies-
FISHworkflow
MicrobeTracker (spotFinderF)

-References-

-Author-
Sangjin Kim, modified 2021, May 12
%}

pixel = 64.2; %nm size of the pixel
load(YFPfilepath);
cellcount = 0; spot = [];
for frame=1:length(spotList)
    for cell=1:length(spotList{frame})
        if ~isempty(spotList{frame}{cell}) 
            cellcount = cellcount + 1;
            spot = [spot; spotList{frame}{cell}.h, spotList{frame}{cell}.w, spotList{frame}{cell}.b, spotList{frame}{cell}.x, spotList{frame}{cell}.y, spotList{frame}{cell}.m, 0];
        end;
    end;
end;

% remove close spots in the single-colored images
for i = 1: size(spot,1)
    for j = 1:(size(spot,1)-i)
        dis = sqrt((spot(i,4)-spot(i+j,4))^2+(spot(i,5)-spot(i+j,5))^2);
        if dis < 10 % within 10 pixels
            spot(i,7) = 2; spot(i+j,7) = 2;     
        end;
    end;
end;
cellcount1 = cellcount-length(find(spot(:,7)==2));
spot(find(spot(:,7)==2),:) = [];
spot1 = spot; %list of spots in YFP channel


load(Cy5filepath);
cellcount = 0; spot = [];
for frame=1:length(spotList)
    for cell=1:length(spotList{frame})
        if ~isempty(spotList{frame}{cell}) 
            cellcount = cellcount + 1;
            spot = [spot; spotList{frame}{cell}.h, spotList{frame}{cell}.w, spotList{frame}{cell}.b, spotList{frame}{cell}.x, spotList{frame}{cell}.y, spotList{frame}{cell}.m, 0];
        end;
    end;
end;
% remove close spots in the single-colored images
for i = 1: size(spot,1)
    for j = 1:(size(spot,1)-i)
        dis = sqrt((spot(i,4)-spot(i+j,4))^2+(spot(i,5)-spot(i+j,5))^2);
        if dis < 10 % within 10 pixels
            spot(i,7) = 2; spot(i+j,7) = 2;     
        end;
    end;
end;
cellcount2 = cellcount-length(find(spot(:,7)==2));
spot(find(spot(:,7)==2),:) = [];
spot2 = spot; %list of spots in Cy5 channel


load(Cy3filepath);
cellcount = 0; spot = [];
for frame=1:length(spotList)
    for cell=1:length(spotList{frame})
        if ~isempty(spotList{frame}{cell}) 
            cellcount = cellcount + 1;
            spot = [spot; spotList{frame}{cell}.h, spotList{frame}{cell}.w, spotList{frame}{cell}.b, spotList{frame}{cell}.x, spotList{frame}{cell}.y, spotList{frame}{cell}.m, 0];
        end;
    end;
end;
% remove close spots in the single-colored images
for i = 1: size(spot,1)
    for j = 1:(size(spot,1)-i)
        dis = sqrt((spot(i,4)-spot(i+j,4))^2+(spot(i,5)-spot(i+j,5))^2);
        if dis < 20 % within 10 pixels
            spot(i,7) = 2; spot(i+j,7) = 2;     
        end;
    end;
end;
cellcount3 = cellcount-length(find(spot(:,7)==2));
spot(find(spot(:,7)==2),:) = [];
spot3 = spot; %list of spots in Cy3 channel

% check spot intensity histogram
% figure, hist(spot1(:,1));

% from the histogram, identify the upper limit of useful spots. 
% if images contain a clump of beads (very bright) as well as
% individual beads, set the upperlimit intensity for individual beads to
% filter out bright spots.
upperlimit = 0.2;
upperlimit = max(spot1(:,1)); %made this line for public distribution.
threshold = 3; % threshold for identifying the same particle 3*64 = 192 nm??

spotmatched = []; L12 = []; L13 = [];
for i = 1:cellcount1
    clear tmp*;
    if spot1(i,1)<upperlimit
        %distance between spot centroid in YFP and Cy5 channels
        spot2tmp = sqrt((spot2(:,4)-spot1(i,4)).^2+(spot2(:,5)-spot1(i,5)).^2);
        %distance between spot centroid in YFP and Cy3 channels
        spot3tmp = sqrt((spot3(:,4)-spot1(i,4)).^2+(spot3(:,5)-spot1(i,5)).^2);
        [tmpL2, tmpi2] = min(spot2tmp);
        [tmpL3, tmpi3] = min(spot3tmp);
        
        if tmpL2 < threshold && tmpL3 < threshold && spot2(tmpi2,7) == 0 && spot3(tmpi3,7) == 0
            spot1(i,7) = 1; spot2(tmpi2,7) = 1; spot3(tmpi3,7) = 1; % flag that it is counted
            %spots that appeared in YFP, Cy5, and Cy3 channels
            spotmatched = [spotmatched; spot1(i,:),spot2(tmpi2,:), spot3(tmpi3,:)]; 
            L12 = [L12; sqrt((spot2(tmpi2,4)-spot1(i,4))^2+(spot2(tmpi2,5)-spot1(i,5))^2)*pixel];
            L13 = [L13; sqrt((spot3(tmpi3,4)-spot1(i,4))^2+(spot3(tmpi3,5)-spot1(i,5))^2)*pixel];
        end;
    end;
end;

%Calculate Affine transformation matrix for YFP-Cy5
%Cy5 channel coordinates 
X = [spotmatched(:,11), spotmatched(:,12), ones(size(spotmatched,1),1)];
%YFP channel coordinates
Y1 = spotmatched(:,4); Y2 = spotmatched(:,5);

%linear regression
A = mvregress(X,Y1);
B = mvregress(X,Y2);

C=A/sqrt(A(1)^2+A(2)^2);
D=B/sqrt(B(1)^2+B(2)^2);
E=(sqrt(A(1)^2+A(2)^2)+sqrt(B(1)^2+B(2)^2))/2;

%calculate residuals
Z1 = X*A-Y1; Z2 = X*B-Y2; Z3 = X*C*E-Y1; Z4 = X*D*E-Y2;
sqrt(A(3)^2+B(3)^2)*pixel % how much shift? translation factor only
Err1 = sqrt( Z1.^2 + Z2.^2)*pixel; 
Err2 = sqrt( Z3.^2 + Z4.^2)*pixel;

figure, hist(Z1,100); figure, hist(Z4, 100);
figure, hist(Err1,100);
% figure, hist(Err2,100);
length(find(Err1<20))/length(Err1) % percentage of the remaining error less than 20 nm
mean(Err1) % remaining error (in nm). 
Xans = A; %ABaffine1
Xans(:,2) = B; %ABaffine2

%Calculate Affine transformation matrix for YFP-Cy3
%Cy3 channel coordinates 
X = [spotmatched(:,18), spotmatched(:,19), ones(size(spotmatched,1),1)];

%linear regression
A = mvregress(X,Y1);
B = mvregress(X,Y2);

C=A/sqrt(A(1)^2+A(2)^2);
D=B/sqrt(B(1)^2+B(2)^2);
E=(sqrt(A(1)^2+A(2)^2)+sqrt(B(1)^2+B(2)^2))/2;

%calculate residuals
Z1 = X*A-Y1; Z2 = X*B-Y2; Z3 = X*C*E-Y1; Z4 = X*D*E-Y2;
sqrt(A(3)^2+B(3)^2)*pixel % how much shift? translation factor only
Err1 = sqrt( Z1.^2 + Z2.^2)*pixel; 
Err2 = sqrt( Z3.^2 + Z4.^2)*pixel;

figure, hist(Z1,100); figure, hist(Z4, 100);
figure, hist(Err1,100);
% figure, hist(Err2,100);
length(find(Err1<20))/length(Err1) % percentage of the remaining error less than 20 nm
mean(Err1) % remaining error (in nm). 

Xans(:,3) = A; %ACaffine1
Xans(:,4) = B; %ACaffine2



