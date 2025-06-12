function optimal_paramspace_app(D, alpha)
    % plot optimal parameter space for variable static localization error
    % optimal_param - Plots colormaps of the relative error of the MSD
    % with variable localization error
    % Inputs:
    %   D        - Diffusion constant (um^2/s^\alpha)
    % alpha      - Anomalous exponent
   
    figure;
    hold on;
    
    % imaging parameter %% CHANGE THIS as needed
    psf = 0.256; %um
    p = 164/0.02; % 164 photons from a single mEos3.2 from tE = 20ms
    
    n=1:1:10; % frame lag
    t = 0.01:0.01:0.2; % frame illumination time (tE range; unit: s)
    
    % gridpoints based on frame lag and illumination time
    [N, T] = meshgrid(n, t);
    result = NaN(size(n));
    
    true = log(2*D*(T.*N).^(alpha)); % original MSD
    
    % MSD with motion averaging and static localization error due to the particle motion
    % Define new S, which stands for sigma^2.
    S = (psf^2) ./ (p .* T) + (2 * D .* T.^alpha) ./ ((alpha + 2) * (alpha + 1) * p .* T);
    msd_mb = log((2*D*T.^(alpha)).*((N+1).^(alpha+2)+(N-1).^(alpha+2)-2*(N).^(alpha+2))/((alpha+2)*(alpha+1))- 4*D*T.^(alpha)/((alpha+2)*(alpha+1))+2*S) ; 
    result= abs(1-(1./true).*msd_mb); %relative error
    
%     result= abs(msd_mb-true); % error in terms of log(MSD)-log(MSD')

    figure;
    imagesc(n, t, result);
    set(gca, 'YDir', 'normal', 'FontSize', 12, 'FontWeight', 'bold');
    colormap("parula");
    colorbar;
    xlabel('Frame lag','FontSize', 14, 'FontWeight', 'bold');
    ylabel('t_{E} (s)', 'FontSize', 14, 'FontWeight', 'bold');
    threshold = 0.05;
    [idxRows, idxCols] = find(abs(result) <= threshold); % relative error would always be less than or equal to 5%
    nPoints = n(idxCols); 
    tPoints = t(idxRows); 
    hold on;
    scatter(nPoints, tPoints, 15, 'white', 'filled'); hold on;
    set(gca, 'FontSize',16,'FontWeight','bold','LineWidth',0.5);hold on;
    yticks([0.02 0.05 0.08 0.11 0.14 0.17 0.2])
end