%% FISH image analysis procedure is explained here step-by-step.
%% DO NOT execute this script entirely. Make sure you read through and use parts one at a time.

%%
% ----------
% 1. Obtain cell meshes by running MicrobeTracker or Oufti
% In the microbeTracker, load both signal1 and signal2 images and perform background subtraction (radio button: all)
% For microbeTracker, use Matlab2013 or Matlab2014. For Oufti, use recent Matlab (e.g., 2019).
% Example parameter set for phase contrast images: Ecoli-largepix.set (160 nm pix). Ecoli-smallpix.set (64.5 nm pix)

%%
% ----------
% 2. Identify spots and obtain their info (e.g., magnitude, from gaussian fitting) by using spotFinderZ in microbeTracker
% spots in channel 1 (Cy5, Z5) and spots in channel 2 (Cy3, Z3) are saved in different cellList mat files.
% Example spot parameter set: Cy5spotparam.sfp and Cy3spotparam.sfp
% To view spots, use dispcellall function in microbeTracker.
images = loadimageseries('folder path to fluorescence images'); 
load('mesh-spot.mat');
dispcellall(cellList,images,[],[],[],[],3,'circle'); 

%%
% ----------
% 3. Obtain threshold for spot curation
SegHist(cellList_timezero); % At time zero, segment intensities indicate background level. 
% From the histogram plot, we take the end of the Gaussian curve, i.e., maximum of the background level as the threshold for real spots.

%%
% ----------
% 4. Find SM normalization factor
% This function makes a list of spots in t1 (or t12 after repression experiment) images, and hence a list of single-molecule spots.
ListSpotmag(folderPath,channel,timePoints,thresh1,thresh2)
% Save the result (Xans, nx1 array) in Excel. Have results from different
% days as independent columns. Top of the column goes the description of the column e.g.,
% "day1z5spot" "day1z3spot". 
% Run GMM_Fitting in python to run gmm (gaussian mixture model)
% system('python GMM_Fitting.py');

%%
% ----------
% 5. smFISH based on spot intensity only
% mRNA number per cell is calculated from the sum of spot intensities
% within a cell
masterfolder = 'J:\Shares\Data_02\Sangjin Kim\Microscope Data\FISH trials\12early_lacZ\';
dayList(1).folderpath = strcat(masterfolder,'140505-lacZ 2CFISH IPTG glucose pulse burstcheck-NSTORM\B.MGhighIPTGearly\');
dayList(2).folderpath = strcat(masterfolder,'140611-lacZ 2CFISH IPTG glucose pulse-NSTORM\B.MGhighearly\');
dayList(3).folderpath = strcat(masterfolder,'140613-lacZ 2CFISH IPTG glucose pulse high-NSTORM\B.MGhighearly\');
dayList(4).folderpath = strcat(masterfolder,'140716-lacZ 2CFISH IPTG glucose pulse high-NSTORM\A.MGhighearly\');

% threshold for spot curation (step #3)
daily_thresh1(1) = 0.003; daily_thresh2(1) = 0.005; 
daily_thresh1(2) = 0.0025; daily_thresh2(2) = 0.007;
daily_thresh1(3) = 0.003; daily_thresh2(3) = 0.006;
daily_thresh1(4) = 0.003; daily_thresh2(4) = 0.006;
daily_thresh = [daily_thresh1', daily_thresh2'];

% SM normalization factor from above process (step #4)
daily_scaleZ5(1) = 0.36494524; daily_scaleZ3(1) = 0.49707307;
daily_scaleZ5(2) = 0.28242603; daily_scaleZ3(2) = 0.55073717;
daily_scaleZ5(3) = 0.29541658; daily_scaleZ3(3) = 0.49020809;
daily_scaleZ5(4) = 0.34345343; daily_scaleZ3(4) = 0.67065377;
daily_scale = [daily_scaleZ5', daily_scaleZ3'];

% time course information
% time (minutes) of sampling (FISH )
timecourseInfo.time = [0,1,2,3,4,5,6,7,8,9,10,12]; 
% days for each time point. For example t = 12 is taken in day 2, 3, 4
timecourseInfo.day = {[1,2,3,4],[1,2,3,4],[1,2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[3,4],[2,3],[2,3,4]};
% file name corresponding to the timepoint in each day. For example, t = 12
% is taken in day2-t12, day3-t12, day4-t12.
timecourseInfo.file = {[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4],[5,5,5],[6,6,6],[7,7,7],[8,8,8],[9,9,9],[10,10],[11,11],[12,12,12]}; 

Xans = RNApercell_by_spotmag_bootstrap(dayList,daily_scale,daily_thresh,timecourseInfo);
%Xans(:,1:3) mean, std, fano factor of Z5 per cell
%Xans(:,4:6) mean, std, fano factor of Z3 per cell
%Xans(:,7) total cell number analyzed
%*** parfor loop can create error. Change "parfor" to "for", but bootstrap
%with "for" loop can take a long time.

%%
% ----------
% 6 (opt). smFISH based on total fluorescence signal
% mRNA number per cell is calculated from the total fluorescence intensity
% within the cell boundary, (subtracted with the basal intensity level)
% normalized by single-molecule fluorescence intensity.

masterfolder = 'J:\Shares\Data_02\Sangjin Kim\Microscope Data\FISH trials\12early_lacZ\';
dayList(1).folderpath = strcat(masterfolder,'140505-lacZ 2CFISH IPTG glucose pulse burstcheck-NSTORM\B.MGhighIPTGearly\');
dayList(2).folderpath = strcat(masterfolder,'140611-lacZ 2CFISH IPTG glucose pulse-NSTORM\B.MGhighearly\');
dayList(3).folderpath = strcat(masterfolder,'140613-lacZ 2CFISH IPTG glucose pulse high-NSTORM\B.MGhighearly\');
dayList(4).folderpath = strcat(masterfolder,'140716-lacZ 2CFISH IPTG glucose pulse high-NSTORM\A.MGhighearly\');

% SM normalization factor from above process (step #4)
daily_scaleZ5(1) = 0.36494524; daily_scaleZ3(1) = 0.49707307;
daily_scaleZ5(2) = 0.28242603; daily_scaleZ3(2) = 0.55073717;
daily_scaleZ5(3) = 0.29541658; daily_scaleZ3(3) = 0.49020809;
daily_scaleZ5(4) = 0.34345343; daily_scaleZ3(4) = 0.67065377;
daily_scale = [daily_scaleZ5', daily_scaleZ3'];

% Get average fluorescence pixel value at time zero = background fluroescence
% per pixel by meaninthist function in microbeTracker.
%e.g., x(:,1) = meaninthist(cellList,'signal1'); x(:,2) = meaninthist(cellList,'signal2');
% Next, use gmm code (in step #4) GMM_Fitting.py with k = 1, one gaussian fit
% to obtain average base pixel value for Cy5 and Cy3 channel.
daily_avgbasePix1(1) = 0.0005091; daily_avgbasePix2(1) = 0.00071035;
daily_avgbasePix1(2) = 0.0005255; daily_avgbasePix2(2) = 0.00079067;
daily_avgbasePix1(3) = 0.00061279; daily_avgbasePix2(3) = 0.00116339;
daily_avgbasePix1(4) = 0.00049022; daily_avgbasePix2(4) = 0.00062315;
daily_avgbasePix = [daily_avgbasePix1', daily_avgbasePix2'];

% time course information (same as step 5)
% time (minutes) of FISH 
timecourseInfo.time = [0,1,2,3,4,5,6,7,8,9,10,12];
% days for each time point
timecourseInfo.day = {[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[3,4],[2,3],[2,3,4]};
% file name corresponding to the timepoint in each day
timecourseInfo.file = {[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4],[5,5,5],[6,6,6],[7,7,7],[8,8,8],[9,9,9],[10,10],[11,11],[12,12,12]};

%In the following algorithm, sum of fluorescence values within the cell is
%subtracted with cell area times avgbasePixel value and then divided by SM
%normalization factor.
Xans = RNApercell_by_totalfluorescence_bootstrap(dayList,daily_scale,daily_avgbasePix,timecourseInfo);
%Xans(:,1:3) mean, std, fano factor of Z5 per cell
%Xans(:,4:6) mean, std, fano factor of Z3 per cell
%Xans(:,7) total cell number analyzed

%%
% ----------
% 7. To estimate premature termination
% example data: Promoter turn-off at t = 90 s after 0.2 mM IPTG induction
% (one day example)
% mRNA lifetime (mean) was 1.52 min (Z5) and 1.66 min (Z3)
X = [0,0.0134487000000000,0.0181838000000000;1,1.79951000000000,0.0261977000000000;2,4.73399000000000,0.122877000000000;3,3.79987000000000,0.287993000000000;4,3.08398000000000,1.61866000000000;5,1.32696000000000,1.80267000000000;6,0.976509000000000,1.54928000000000;7,0.587405000000000,0.963588000000000;8,0.288642000000000,0.430795000000000;9,0.206146000000000,0.236055000000000;10,0.172012000000000,0.194890000000000;12,0.104088000000000,0.0946338000000000];
intP = IntegrationFISH(X,[1.52,1.66]);

%%
% ----------
% 8. Co-localization analysis
% first, get Affine transformation matrix from tetraspek bead images
% bead samples are imaged in YFP, Cy5, and Cy3 channels (no cells) using
% the same microscope. 
% Bead images are analyzed by spotFinderF of microbeTracker, which yield "spotList" file containing spot locations and
% intensities.
% Make sure to verify the spots by overlaying with their fluorescence images.
images = loadimageseries('folder path to fluorescence images'); 
load('spot file.mat');
dispspotsall(spotList,images,'circle'); 

% filepath = spotList file.
Affine_matrix = getAffineFactor(YFPfilepath, Cy5filepath, Cy3filepath);
% Get Affine_matrix for each day.

%----
masterfolder = 'J:\Shares\Data_02\Sangjin Kim\Microscope Data\FISH trials\12early_lacZ\';
dayList(1).folderpath = strcat(masterfolder,'140505-lacZ 2CFISH IPTG glucose pulse burstcheck-NSTORM\B.MGhighIPTGearly\');
dayList(2).folderpath = strcat(masterfolder,'140611-lacZ 2CFISH IPTG glucose pulse-NSTORM\B.MGhighearly\');
dayList(3).folderpath = strcat(masterfolder,'140613-lacZ 2CFISH IPTG glucose pulse high-NSTORM\B.MGhighearly\');
dayList(4).folderpath = strcat(masterfolder,'140716-lacZ 2CFISH IPTG glucose pulse high-NSTORM\A.MGhighearly\');

% threshold for spot curation (step #3)
daily_thresh1(1) = 0.003; daily_thresh2(1) = 0.005; 
daily_thresh1(2) = 0.0025; daily_thresh2(2) = 0.007;
daily_thresh1(3) = 0.003; daily_thresh2(3) = 0.006;
daily_thresh1(4) = 0.003; daily_thresh2(4) = 0.006;
daily_thresh = [daily_thresh1', daily_thresh2'];

% SM normalization factor from above process (step #4)
daily_scaleZ5(1) = 0.36494524; daily_scaleZ3(1) = 0.49707307;
daily_scaleZ5(2) = 0.28242603; daily_scaleZ3(2) = 0.55073717;
daily_scaleZ5(3) = 0.29541658; daily_scaleZ3(3) = 0.49020809;
daily_scaleZ5(4) = 0.34345343; daily_scaleZ3(4) = 0.67065377;
daily_scale = [daily_scaleZ5', daily_scaleZ3'];

% daily_affine from Affine_matrix (each columns are ABaffine1, ABaffine2,
% ACaffine1, ACaffine2)
daily_affine(1).ABaffine1 =[1.00258547845869;-0.000145160554684541;-1.60658141374572;];
daily_affine(1).ABaffine2 =[2.17330179733176e-05;1.00283409924542;-2.09966855254672;];
daily_affine(1).ACaffine1 =[1.00086354655070;1.63014702434462e-05;-0.626201955735291;];
daily_affine(1).ACaffine2 =[0.000140414796535213;1.00081230890480;-0.862637089350683;];
            
daily_affine(2).ABaffine1 = [1.00337046063806;-1.01676478540238e-05;-2.25909239749171;];
daily_affine(2).ABaffine2 = [-4.52539199503772e-05;1.00333242660941;-2.37287366564344;];
daily_affine(2).ACaffine1 = [1.00130639753223;6.55294641219577e-05;-0.966481742065977;];
daily_affine(2).ACaffine2 = [0.000129691632672423;1.00129046496446;-1.02786580354267;];

daily_affine(3).ABaffine1 = [1.00319374608803;-4.20221763038701e-06;-2.07992686189881;];
daily_affine(3).ABaffine2 = [4.23964591229039e-05;1.00330543205515;-2.45836721797975;];
daily_affine(3).ACaffine1 = [1.00112078302782;8.85685196178556e-05;-0.881879461965597;];
daily_affine(3).ACaffine2 = [5.64992947825514e-05;1.00105527241642;-0.855695595934417;];

daily_affine(4).ABaffine1 = [1.00342879273967;-6.25902825057270e-05;-2.38441412518357;];
daily_affine(4).ABaffine2 = [-2.71378010821591e-05;1.00355662772723;-2.13404951627375;];
daily_affine(4).ACaffine1 = [1.00125805865084;6.92988031554060e-05;-0.989591800213048;];
daily_affine(4).ACaffine2 = [0.000124515243279087;1.00117053119993;-0.803622158255415;];

% time course information
% time (minutes) of sampling (FISH )
timecourseInfo.time = [0,1,2,3,4,5,6,7,8,9,10,12]; 
% days for each time point. For example t = 12 is taken in day 2, 3, 4
timecourseInfo.day = {[1,2,3,4],[1,2,3,4],[1,2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[3,4],[2,3],[2,3,4]};
% file name corresponding to the timepoint in each day. For example, t = 12
% is taken in day2-t12, day3-t12, day4-t12.
timecourseInfo.file = {[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4],[5,5,5],[6,6,6],[7,7,7],[8,8,8],[9,9,9],[10,10],[11,11],[12,12,12]}; 

Xans = RNAspot_bootstrap(dayList,daily_scale,daily_thresh,daily_affine,timecourseInfo);
%parfor loop can create error. Chage parfor to for, but bootstrap with for
%can take a long time.

