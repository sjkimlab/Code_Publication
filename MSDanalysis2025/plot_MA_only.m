function plot_MA_only(X0, Y0, X1, Y01, t, truetraj,  trajMA) 
%Figure1b
    figure,
    lag_times_avg = exp(X1);
    plot(X0, Y0, 'LineWidth', 4,'LineStyle','-','Marker','square','color','[0.5, 0.5, 0.5]');hold on; % true MSD
    plot(X1, Y01, 'LineWidth', 4,'LineStyle','-','Marker','o','color','[1, 0, 0]');hold on; % MA only MSD
    xlabel('log(lag time (s))');
    ylabel('log(MSD (\mum^{2}))');
    set(gca,'linewidth',2, 'fontsize', 12);
    legend('True', 'MA');
    pbaspect([1 1 1]);
    axes('Position',[.7 .7 .2 .2]);
    hold on;
    box on
    plot(t, truetraj, 'LineWidth', 2, 'color','[0.5, 0.5, 0.5]'); hold on; % true trajectory
    plot(lag_times_avg,  trajMA, 'LineWidth', 2, 'color','[1, 0, 0]');hold on; %  MA only trajectory
    xlabel('Time (s)');
    ylabel('Position (\mum^2)');
    set(gca,'linewidth',2, 'fontsize', 12);
end