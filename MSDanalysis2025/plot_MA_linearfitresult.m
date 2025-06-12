function plot_MA_linearfitresult(Y1,Y2, X1, X2, ...
                           D, alpha)
    %    D            - initial guess diffusion coefficient (um^2/s^alpha)
    % alpha           - initial guess diffusion exponent 
    %   Y1            - log of MSD'(MA1+LocE)
    %   Y2            - log of MSD'(MA2+LocE)
    %   X1            - log of lag times (MA1)
    %   X2            - log of lag times (MA2)

    nl1 = floor(length(X1)/4); % We will take first 25% of the datapoints
    nl2 = floor(length(X2)/4);
    msd_lag_MB1_avg = exp(Y1); % MSD' values for MA1
    msd_lag_MB2_avg = exp(Y2); % MSD' values for MA2
    lag_times_avg = exp(X1);   % lag time values for MA1
    lag_times_avg2 = exp(X2);  % lag time values for MA2
    
    % Extract linear model fitted parameters
    [alpha1, D1] = estimate_diffusion_constant_and_exponent(msd_lag_MB1_avg, lag_times_avg);
    legendText1 = sprintf('$$D_{e} = %.2f,\\ \\alpha_{e} = %.2f$$', D1, alpha1);

    [alpha2, D2] = estimate_diffusion_constant_and_exponent(msd_lag_MB2_avg, lag_times_avg2);
    legendText2 = sprintf('$$D_{e} = %.2f,\\ \\alpha_{e} = %.2f$$', D2, alpha2);

    % Plotting
    fig = figure;
    hold on;

    plot(X1, Y1, 'LineWidth', 2, 'Marker', 'o', 'LineStyle', 'none');
    plot(X1(1:nl1), log(2 * D1 .* lag_times_avg(1:nl1).^alpha1), '--', 'LineWidth', 2);
    plot(X2, Y2, 'LineWidth', 2, 'Marker', 'square', 'LineStyle', 'none'); 
    plot(X2(1:nl2), log(2 * D2 .* lag_times_avg2(1:nl2).^alpha2), '--', 'LineWidth', 2);
    box on;

    xlabel('log(lag time (s))', 'Interpreter', 'latex');
    ylabel('log(MSD ($\mu m^{2}$))', 'Interpreter', 'latex');
    set(gca, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 2);
    legend({ 'MSD$^\prime$ ($\gamma = t_E$)', legendText1, ...
             'MSD$^\prime$ ($\gamma = 5t_E$)', legendText2}, ...
            'Interpreter', 'latex', 'FontSize', 12, 'Location', 'southeast');
    legend boxoff;
    pbaspect([1 1 1]);
end