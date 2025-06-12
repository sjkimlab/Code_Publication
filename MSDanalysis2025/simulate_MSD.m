function [Y0, Y1, Y01, Y2, truetraj, trajMA, trajMALocE] = simulate_MSD(T, dt, D, alpha, locError, num_trajectories, Tstop, frame_time)
    % simulate_MSD - Generate FBM trajectories and compute EA-TA-MSD curves
    % 
    % Inputs:
    %   T                 - total observation time (s)
    %   dt                - simulation time step (s)
    %   D                 - diffusion coefficient (um^2/s^alpha)
    %   alpha             - anomalous exponent
    %   locError          - localization error (um)
    %   num_trajectories  - number of trajectories
    %   Tstop             - camera off time for MA2 method (s)
    %   frame_time        - frame acquisition time (s)
    %
    % Outputs:
    %   Y0            - log of True MSD (um2)
    %   Y1            - log of MSD'(MA1+LocE) for Tstop = 0
    %   Y01           - log of MSD' (only MA)
    %   Y2            - log of MSD'(MA2+LocE) based on input Tstop
    % truetraj        - true trajectory
    % trajMA          - trajectory with MA1
    % trajMALocE      - trajectory with both MA1 and LocE

    num_steps = round(T/dt);   
    t = dt*(1:num_steps);   

    % Construct covariance matrix
    covMat = zeros(num_steps);
    for i = 1:num_steps
        for j = 1:i
            covMat(i, j) = t(i)^alpha + t(j)^alpha - abs(t(i)-t(j))^alpha;
        end
    end

    % Cholesky decomposition
    L = chol(D * covMat, 'lower');

    % Generate FBM trajectories
    tf(num_trajectories, 1).traj = [];
    for i = 1:num_trajectories
        tf(i).traj = L * randn(num_steps, 1);
    end
    truetraj = tf(1).traj;
    max_lag = num_steps - 1; 
    lag_times = (1:max_lag) * dt;
    Tstop_interval = floor(Tstop/dt);
    interval = floor(frame_time / dt);
    num_intervals = floor(T / frame_time);
    num_intervals2 = floor(T / (frame_time + Tstop));
    lag_times_avg = (1:num_intervals - 1) * frame_time;
    lag_times_avg2 = (Tstop + frame_time) + ((1:num_intervals2 - 1) - 1) * (Tstop + frame_time);

    locE = normrnd(0, locError, [num_intervals, 1]);
    locE1 = locE(1:num_intervals2);

    averaged_trajectory = zeros(num_trajectories, num_intervals);
    trajectory_MB1 = zeros(num_trajectories, num_intervals);
    averaged_traj_MB2 = zeros(num_trajectories, num_intervals2);
    trajectory_MB2 = zeros(num_trajectories, num_intervals2);
    msd_lag_all = zeros(num_trajectories, max_lag);
    msd_lag_MB1_all = zeros(num_trajectories, num_intervals - 1);
    msd_lag_MB2_all = zeros(num_trajectories, num_intervals2 - 1);
    msd_lag_mb_all = zeros(num_trajectories, num_intervals - 1);

   for traj_idx = 1:num_trajectories
        trajectory = tf(traj_idx).traj;

        % True MSD calculation
        for lag = 1:max_lag           
            disp = trajectory(lag+1:end) - trajectory(1:end-lag);
            msd_lag_all(traj_idx, lag) = mean(disp.^2);
        end

        % MA1- Motion-avergaing effect - contineous exposure
        for i = 1:num_intervals
            idx_range = (i - 1) * interval + (1:interval);
            averaged_trajectory(traj_idx, i) = mean(trajectory(idx_range));      % MA1 
            trajectory_MB1(traj_idx, i) = averaged_trajectory(traj_idx, i) + locE(i); % MA1 + LocE 
        end

        for lag_avg = 1:num_intervals - 1     
            displacements_mb = averaged_trajectory(traj_idx, lag_avg + 1:end) - averaged_trajectory(traj_idx,1:end - lag_avg);
            displacements_MB1 = trajectory_MB1(traj_idx, lag_avg + 1:end) - trajectory_MB1(traj_idx, 1:end - lag_avg);
            msd_lag_MB1_all(traj_idx, lag_avg) = mean(displacements_MB1.^2); % TA-MSD (MA1)
            msd_lag_mb_all(traj_idx, lag_avg) = mean(displacements_mb.^2);   % TA-MSD MA only
        end

        % MA2 - Motion-avergaing effect - time lapse exposure
        for i = 1:num_intervals2
            idx_range1 = (i - 1) * (interval + Tstop_interval) + (1:interval); 
            % Time averaging (TA)
            averaged_traj_MB2(traj_idx, i) = mean(trajectory(idx_range1));     % MA2
            trajectory_MB2(traj_idx, i) = averaged_traj_MB2(traj_idx, i) + locE1(i);  % MA2 + LocE         
        end

        for lag_avg2 = 1:num_intervals2 - 1    
            % Time averaging (TA)
            displacements_MB2 = trajectory_MB2(traj_idx, lag_avg2 + 1:end) - trajectory_MB2(traj_idx, 1:end - lag_avg2); 
            msd_lag_MB2_all(traj_idx, lag_avg2) = mean(displacements_MB2.^2);  % TA-MSD (MA2)
        end
   end

   trajMA = averaged_trajectory(1,1:end-1);
   trajMALocE = trajectory_MB1(1,1:end-1);

    % Ensemble averaging (EA)
    msd_lag_avg = mean(msd_lag_all, 1);             % True MSD
    msd_lag_mb_avg = mean(msd_lag_mb_all, 1);       % MA only
    msd_lag_MB1_avg = mean(msd_lag_MB1_all, 1);     % MA1+locE
    msd_lag_MB2_avg = mean(msd_lag_MB2_all, 1);     % MA2+locE 

    % log(MSD) and log(MSD')
    Y0 = log(msd_lag_avg);
    Y01 = log(msd_lag_mb_avg);
    Y1 = log(msd_lag_MB1_avg);
    Y2 = log(msd_lag_MB2_avg);
end