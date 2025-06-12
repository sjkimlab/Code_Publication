function plot_MA_locE(X0, Y0, X1, Y01, Y1, trajMA, trajMALocE)
% Figure1c
    figure,
    lag_times_avg = exp(X1);
    plot(X0, Y0, 'LineWidth', 0.5,'LineStyle','-');hold on; % true MSD
    plot(X1, Y01, 'LineWidth', 2,'Marker','o','color','[1, 0, 0]');hold on; % MSD' (MA only MSD
    plot(X1, Y1, 'LineWidth', 2,'Marker','diamond','color','[0, 0.45, 0.74]');hold on;  % MSD' (MA+LocE) 
    xlabel('log(lag time (s))'); 
    ylabel('log(MSD (\mum^{2}))');
    set(gca,'linewidth',2, 'fontsize', 12);
    legend('True','MA', 'MA+LocE');
    pbaspect([1 1 1]);
    axes('Position',[.7 .7 .2 .2]);
    hold on;
    box on
    plot(lag_times_avg, trajMA, 'LineWidth', 2);hold on; % MB only trajectory
    plot(lag_times_avg, trajMALocE, 'LineWidth', 2);hold on;       % MB+LocE trajectory
    xlabel('Time (s)');
    ylabel('Position (\mum^2)');
    set(gca,'linewidth',2, 'fontsize', 12);
end