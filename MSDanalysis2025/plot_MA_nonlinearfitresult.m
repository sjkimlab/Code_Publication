function plot_MA_nonlinearfitresult(Y1,Y2, X1, X2, ...
                           D, alpha)

    % input parameters
    %   Y1            - log of MSD'(MA1+LocE) for Tstop = 0
    %   Y2            - log of MSD'(MA2+LocE) based on input Tstop
    %   X1            - log of lag times for MA1
    %   X2            - log of lag times for MA2    
    %   D, alpha      - initial guess

    nl1 = floor(length(X1)/4); %fit first 25%
    nl2 = floor(length(X2)/4);
    
     % Fit using fsimple model
    [bestparams1, fitted_f1] = fsimple(X1(1:nl1), Y1(1:nl1), D, alpha);
    [bestparams2, fitted_f2] = fsimple(X2(1:nl2), Y2(1:nl2), D, alpha);

    % Extract fitted parameters from nonlinear fitting method
    De1 = bestparams1.a;
    alphae1 = bestparams1.b;
    legendText1 = sprintf('$$D_{e} = %.2f,\\ \\alpha_{e} = %.2f$$', De1, alphae1);

    De2 = bestparams2.a;
    alphae2 = bestparams2.b;
    legendText2 = sprintf('$$D_{e} = %.2f,\\ \\alpha_{e} = %.2f$$', De2, alphae2);

    % Plotting
    fig = figure;
    hold on;

    plot(X1, Y1, 'LineWidth', 2, 'Marker', 'o', 'LineStyle', 'none');
    plot(X1(1:nl1), fitted_f1, '-', 'LineWidth', 2);
    plot(X2, Y2, 'LineWidth', 2, 'Marker', 'square', 'LineStyle', 'none');
    plot(X2(1:nl2), fitted_f2, '-', 'LineWidth', 2);

    xlabel('log(lag time (s))', 'Interpreter', 'latex');
    ylabel('log(MSD ($\mu m^{2}$))', 'Interpreter', 'latex');
    set(gca, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 2);
    legend({ 'MSD$^\prime$ ($\gamma = t_E$)', legendText1, ...
             'MSD$^\prime$ ($\gamma = 5t_E$)', legendText2}, ...
            'Interpreter', 'latex', 'FontSize', 12, 'Location', 'southeast');
    legend boxoff;
    pbaspect([1 1 1]);
end