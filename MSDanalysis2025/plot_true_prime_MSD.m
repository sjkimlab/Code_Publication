function plot_true_prime_MSD(D, alpha,frame_time, locError)
    % Figure 3(a-d)
    % theoretical MSD' with true MSD
   
    dt = 0.001; % sampling time (s)
    T = 10; % total observation time (s)
    %frame_time  = 0.020; % camera acquisition time(s), equivalent to tE
    num_steps = T/dt;   
    t = dt*(1:num_steps);  
    max_lag = num_steps-1; 
    lag_times = (1:max_lag) * dt;
    num_intervals = floor(T / frame_time);

    array = 1:num_intervals-1; 

    lag_times = (1:max_lag) * dt;
    lag_times_avg = (1:num_intervals - 1) * frame_time; % lag times for MA1 method (gamma = tE)

    X0 = log(lag_times);
    X1 = log(lag_times_avg);

    lag_times = exp(X0);
    % Theoretical MSD'
    YMB1 = (2*D*frame_time^(alpha))*((array+1).^(alpha+2)+(array-1).^(alpha+2)-2*(array).^(alpha+2))/((alpha+2)*(alpha+1))- 4*D*frame_time.^(alpha)/((alpha+2)*(alpha+1)) +2*locError.^2;
    YMBL1 = log(YMB1); 
    
    figure,
    
    plot(X0, log(2*D.*lag_times.^alpha), 'LineWidth', 2,'LineStyle','-');hold on; % true MSD (MSD)
    plot(X1, YMBL1, 'linewidth', 2,'LineStyle','--');hold on; % predicted MSD (MSD')
    xlim([-4 2]);
    xlabel('log(lag time (s))', 'Interpreter', 'latex');
    ylabel('log(MSD ($\mu m^{2}$))', 'Interpreter', 'latex');
    set(gca,'linewidth',2, 'fontsize', 14);
    legend('True','MA+LocE');
    legend box off;
    pbaspect([1 1 1]);
end