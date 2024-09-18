%This is a script to fit either a sinlge-population Guassian or a
%two-population Gaussian to D histogram data. 

%This script was writtin by Laura Troyer for use in Dr. Sangjin Kim's Lab
%at UIUC. This was last edited on Sept. 3, 2024.

%use these parameters to change features of the histogram/fit:
scale = 0.011; %size of bin Width to use for histogram plotting
numPopulations = 2; %use 1 for 1-population fit, 2 for 2-population fit

 %load in the tracksFinal structure containig the diffusion data
 strain = input("strain number to plot: ", "s");
 if strain(1:2) ~= 'SK'
     if strain(1:2) == 'sk'
         strain(1:2) = 'SK';
     else
         strain = ['SK' strain];
     end
 end
 
 MainFolder = 'C:\Users\ltroyer2\Documents\Lab Data\Track and Cell Variables\'; %Make this the folder where the strain data is stored
 oldFolder = cd(MainFolder); %save the workers previous working directory and make the MainFolder directory the active directroy
 files = dir([MainFolder strain '*_tracksFinal*newAll2.mat']);
 if length(files)> 1
     FileChosen = 0;
     while FileChosen == 0
         for i = 1:length(files)
             FileChosen = input(['Use file (0 or 1): ' files(i).name ' ']);
             if FileChosen == 1
                 fileName = files(i).name;
                 break
             end
         end
     end
 else
     fileName = files.name;
 end
 load([MainFolder fileName], 'tracksFinal');
 
 %now select the diffusion data
 data = [tracksFinal.D]*10^12;



%Now plot the histogram and get the histogram data
figure
hold on

 h = histogram(data, 'BinWidth', scale,  'DisplayName', [ strain ', ' 10 'experimental data' 10 'mean = ' num2str(mean(data))]);
 %get the histogram values for plotting later   
 NC = h.Values;
 binEdges = h.BinEdges;
 bM = movmean(binEdges, 2);
 bM = bM(2:end);
 tC = length(h.Data);

 %insert figure title, legend, and axis labels
 title([strain ', Bin Width = ' num2str(scale)]);   
 grid on
 ax = gca;
 ax.FontSize = 12;
 xlabel('D (\mum^2/s)', 'FontSize', 16)
 ylabel('Counts', 'FontSize', 16)
 legend('Location', 'northeast')
 legend show
 xlim([-0.1, 0.5])


%Now get either the 1-population or 2-population Gaussian fits
if numPopulations == 1 %use 1-population fit
    %set up initial guesses and fit
    paramInitials = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0.0001 -0.001 0.01], 'Upper', [Inf 0.2 Inf], 'Startpoint', [20, 0.05, 0.02]);
    ft = fittype('Gaussian(x, c, mu, sigma)', 'independent', {'x'}, 'coefficients', {'c', 'mu', 'sigma'}, 'Options', paramInitials);
    [coeff, gof] = fit(bM', NC', ft) %Note: this line is not suppressed to read out 95 percent confidence intervals on fitting
    rsquare = gof.rsquare;

    %set up variables to plot fit output on top of the histogram to show fit    
    bM_fit = min(bM):0.0011:max(bM);
    yFit = Gaussian(bM_fit, coeff.c, coeff.mu, coeff.sigma);
    
    %plot fit
    p1 = plot(bM_fit, yFit, '-', 'Color', [0 0 0], 'LineWidth', 4, 'DisplayName', ['Guassian Fit' 10 'r^2 = ' num2str(rsquare) 10 'c = ' num2str(coeff.c) ', \mu = ' num2str(coeff.mu) ', \sigma = ' num2str(coeff.sigma)]);

elseif numPopulations == 2 %use 2-population fit
    %set up initial guesses and fit
    paramInitials = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0.001 0 0 0 0 0.04], 'Upper', [inf 1 1 1 20 20], 'Startpoint', [26.5, 0.7, 0.05, 0.02, 0.1, 1]); %startpoints are in order of listed coefficients in next line
    ft = fittype('Gaussian2(x, c, a, mu, sigma, mu2, sigma2)', 'independent', {'x'}, 'coefficients', {'c', 'a', 'mu', 'sigma', 'mu2', 'sigma2'}, 'Options', paramInitials)
    [coeff, gof] = fit(bM', NC', ft)
    rsquare = gof.rsquare;
    
    %set up variables to plot fit output on top of the histogram to show fit 
    bM_fit = min(bM):0.0011:max(bM);
    yFit = Gaussian2(bM_fit, coeff.c, coeff.a, coeff.mu, coeff.sigma, coeff.mu2, coeff.sigma2);
    yFit1 = Gaussian(bM_fit, coeff.a*coeff.c, coeff.mu, coeff.sigma);
    yFit2 = Gaussian(bM_fit, ((1-coeff.a)*coeff.c), coeff.mu2, coeff.sigma2);
    
    %plot fit
    p1 = plot(bM_fit, yFit, '-', 'Color', [0 0 0], 'LineWidth', 4, 'DisplayName', ['Double guassian Fit' 10 'r^2 = ' num2str(rsquare) 10 'c = ' num2str(coeff.c)]);
    p2 = plot(bM_fit, yFit1, ':', 'LineWidth', 4, 'DisplayName', ['pop_1 = ' num2str(coeff.a) ', \mu = ' num2str(coeff.mu) ', \sigma = ' num2str(coeff.sigma)]);
    p3 = plot(bM_fit, yFit2, ':', 'LineWidth', 4, 'DisplayName', ['pop_2 = ' num2str((1-coeff.a)), ', \mu2 = ' num2str(coeff.mu2) ', \sigma2 = ' num2str(coeff.sigma2)]);
end

cd(oldFolder)



