function plot_MSD(T, X0, X1, X2, Y0, Y1, Y2, D, frame_time, alpha)
    % plot of Figure 2a, Figure 2b, Figure 2d, Figure 2e
    % comparision of true MSD with MSD' 
    % Inputs:
    %   X0, Y0   - lag times and True MSD
    %   X1, Y1   - Estimated MSD' with gamma = t_E
    %   X2, Y2   - Estimated MSD' with gamma = 5t_E
    %   D        - Diffusion constant (um^2/s^\alpha)
    % frame_time - camera acquisition time (s) same as tE
    % alpha      - Anomalous exponent
   
    num_intervals = floor(T / frame_time);

    figure;
    hold on;
    locError = 0.02;  % um %change this if needed

    % Theoretical prediction for MSD' (gamma = t_E)
    array = 1:num_intervals-1; % MA1  
    YMB1 = (2*D*frame_time^(alpha))*((array+1).^(alpha+2)+(array-1).^(alpha+2)-2*(array).^(alpha+2))/((alpha+2)*(alpha+1))- 4*D*frame_time.^(alpha)/((alpha+2)*(alpha+1)) +2*locError.^2;
    YMBL1 = log(YMB1); 

    % Plot simulated MSD 
    plot(X0, Y0, 'LineWidth', 2, 'Marker', 'x', 'LineStyle', 'none'); 
    plot(X1, Y1, 'LineWidth', 2, 'Marker', 'o', 'LineStyle', 'none');
    plot(X2, Y2, 'LineWidth', 2, 'Marker', 'square', 'LineStyle', 'none');
 
    % Plot theoretical curve
    plot(X1, YMBL1, 'LineWidth', 2);
    box on;
    xlabel('log(lag time (s))', 'Interpreter', 'latex');
    ylabel('log(MSD ($\mu m^{2}$))', 'Interpreter', 'latex');
    set(gca, 'LineWidth', 2, 'FontSize', 14);
    % xlim([-3 2]);
 
    legend({'True', ...
            'MSD$^\prime$ ($\gamma = t_E$)', ...
            'MSD$^\prime$ ($\gamma = 5t_E$)', ...
            'MSD$^\prime$ Theory'}, ...
            'Interpreter', 'latex', 'Location', 'southeast');
 
    pbaspect([1 1 1]);
end