function [alpha_0, D_alpha_0] = estimate_diffusion_constant_and_exponent(msd_lag, lag_times)
    a = floor(length(msd_lag)/4); % taking first 25% datapoints
    xvalues = log(lag_times(1:a)); 
    yvalues = log(msd_lag(1:a)); 
    fit = polyfit(xvalues, yvalues, 1);
    alpha_0 = fit(1);
    D_alpha_0 = exp(fit(2)) / 2;   
    % R = fit.Rsquared.Ordinary;
end