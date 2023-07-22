classdef Diffusion_Simulator_Original < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        zsectioningCheckBox          matlab.ui.control.CheckBox
        z_maxEditField               matlab.ui.control.NumericEditField
        z_maxEditFieldLabel          matlab.ui.control.Label
        z_minEditField               matlab.ui.control.NumericEditField
        z_minEditFieldLabel          matlab.ui.control.Label
        oftracksEditField            matlab.ui.control.NumericEditField
        oftracksEditFieldLabel       matlab.ui.control.Label
        ButtonGroup_2                matlab.ui.container.ButtonGroup
        excludepolesButton           matlab.ui.control.RadioButton
        wholecellButton              matlab.ui.control.RadioButton
        localizationerrormEditField  matlab.ui.control.NumericEditField
        localizationerrormEditFieldLabel  matlab.ui.control.Label
        EditField                    matlab.ui.control.NumericEditField
        EditFieldLabel               matlab.ui.control.Label
        Dm2sEditField                matlab.ui.control.NumericEditField
        Dm2sEditFieldLabel           matlab.ui.control.Label
        frametimemsEditField         matlab.ui.control.NumericEditField
        frametimemsEditFieldLabel    matlab.ui.control.Label
        offramesEditField            matlab.ui.control.NumericEditField
        offramesEditFieldLabel       matlab.ui.control.Label
        MotionBlurringCheckBox       matlab.ui.control.CheckBox
        CellLengthmEditField         matlab.ui.control.NumericEditField
        CellLengthmEditFieldLabel    matlab.ui.control.Label
        CellRadiusmEditField         matlab.ui.control.NumericEditField
        CellRadiusmEditFieldLabel    matlab.ui.control.Label
        RunButton                    matlab.ui.control.Button
        ButtonGroup                  matlab.ui.container.ButtonGroup
        MembraneButton               matlab.ui.control.RadioButton
        CytoplasmButton              matlab.ui.control.RadioButton
        TabGroup                     matlab.ui.container.TabGroup
        simulatedcellwithtracksTab   matlab.ui.container.Tab
        Image3                       matlab.ui.control.Image
        MSDvstTab                    matlab.ui.container.Tab
        Image4                       matlab.ui.control.Image
        UpdateButton                 matlab.ui.control.Button
        PlottingOptionsButtonGroup   matlab.ui.container.ButtonGroup
        linearButton                 matlab.ui.control.RadioButton
        loglogButton                 matlab.ui.control.RadioButton
        DhistogramsTab               matlab.ui.container.Tab
        Image5                       matlab.ui.control.Image
        Image2                       matlab.ui.control.Image
        xNormTab                     matlab.ui.container.Tab
        UpdateButton_2               matlab.ui.control.Button
        BinWidthEditFieldLabel       matlab.ui.control.Label
        BinWidthEditField            matlab.ui.control.NumericEditField
        Image                        matlab.ui.control.Image
    end



    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            
            close all %close all currently open figures
            
            %Cytoplasm Simulation
            if app.CytoplasmButton.Value == true
                %setting variables from user inputs
                type = 'Cytoplasm';
                R = app.CellRadiusmEditField.Value;   %cell radius
                L = app.CellLengthmEditField.Value;   %cell length
                N = app.offramesEditField.Value;  %number of frames
                frameTime = app.frametimemsEditField.Value/1000;   %frame time
                nTracks = app.oftracksEditField.Value;         %number of tracks
                loc_err = app.localizationerrormEditField.Value;   %localization error
                D = app.Dm2sEditField.Value;  %diffusion coefficient
                a = app.EditField.Value;  %anomalous diffusion coefficient
                binW = app.BinWidthEditField.Value;    %xNorm graph bin width
                percent = "0";         %percent of particles on the membrane
                motion_blur_timing = floor(frameTime.*1000/3)/1000;    %motion blurring timing
                T_exp = N.*frameTime;     %total time
                twenty_five_percent = floor(N/4);          %25 percent of the number of frames (for D and Alpha calculations)
                if app.MotionBlurringCheckBox.Value == true
                    avg = 1;
                end
                if app.MotionBlurringCheckBox.Value == false          %disabling motion blurring if user chooses
                    avg = 0;
                end
                z_min = app.z_minEditField.Value;
                z_max = app.z_maxEditField.Value;
                parameters = struct('Type',type,'frame_time',frameTime,'motion_blur_microstep', motion_blur_timing,'nFrames',N,'total_time',T_exp,'nTracks',nTracks,'L',L,'R',R,'motion_blur',avg,'D',D,'loc_err',loc_err,'z_min',z_min,'z_max',z_max);
                
                %getting the xNorm plot
                %xNorm_plots_runner('R',R,'L', L, 'N', nTracks,'factor',1.5,'percent',percent,'locErr',loc_err,'binW',binW);
                %xNorm_load = imread('xNorm.png');
                %app.Image.ImageSource = xNorm_load;
                %close
                
                %running cytoplasm simulations (with nTracks number of tracks) and putting the output D-values into a histogram
                n_tracks = zeros(1, nTracks);
                for i = 1:nTracks
                    n_tracks(i) = nTracks;
                end
                
                clear tracksFinal
                
                tracksFinal = struct([]);
                xNorm_values = [];
                cc_param = struct('R',R,'L',L-2.*R);
                plot_model(cc_param);
                figure(1)
                
                for ii = 1:n_tracks
                    %T_exp = nFrames(ii).*frameTime;
                    params = struct('dt',motion_blur_timing,'dt_out', frameTime,'t_fin',T_exp,'l0',L,'w0',2.*R,'totR',1,'avg',avg,'D',D,'loc',1,'sigma',loc_err,'boundary',true);
                    [track,params] = rand_diff_3D_SphCyl(params);
                    track(:,[1 2]) = track(:,[2 1]);
                    tracksFinal(ii).pole = any(abs(track(:,2))>=L/2-R);
                    tracksFinal(ii).z = any(track(:,3)>=z_max) + any(track(:,3)<=z_min);
                    if app.excludepolesButton.Value == true
                        if app.zsectioningCheckBox.Value == true
                            if tracksFinal(ii).pole == 0
                                if tracksFinal(ii).z == 0
                                    tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                                    tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                                    xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                                    figure(1)
                                    %hold on
                                    plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                                    daspect([1 1 1])
                                    set(gcf,'color','w');
                                    set(gca,'linewidth',2)
                                    hold on
                                end
                            end
                        end
                        if app.zsectioningCheckBox.Value == false
                            if tracksFinal(ii).pole == 0
                                tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                                tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                                xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                                figure(1)
                                %hold on
                                plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                                daspect([1 1 1])
                                set(gcf,'color','w');
                                set(gca,'linewidth',2)
                                hold on
                            end
                        end
                    end
                    if app.wholecellButton.Value == true
                        if app.zsectioningCheckBox.Value == true
                            if tracksFinal(ii).z == 0
                                tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                                tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                                xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                                figure(1)
                                %hold on
                                plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                                daspect([1 1 1])
                                set(gcf,'color','w');
                                set(gca,'linewidth',2)
                                hold on
                            end
                        end
                        if app.zsectioningCheckBox.Value == false
                            tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                            tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                            xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                            figure(1)
                            %hold on
                            plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                            daspect([1 1 1])
                            set(gcf,'color','w');
                            set(gca,'linewidth',2)
                            hold on
                        end
                    end
                end
                camup([0 1 0])
                xlabel('x (µm)')
                ylabel('y (µm)')
                zlabel('z (µm)')
                savefig('cytoplasm_tracks.fig');
                exportgraphics(gca,'cytoplasm_tracks.jpg','Resolution',300);
                close
                %displaying simulated tracks
                cytoplasm_tracks_load = imread('cytoplasm_tracks.jpg');
                app.Image3.ImageSource = cytoplasm_tracks_load;
                delete('cytoplasm_tracks.jpg')
                
                if app.excludepolesButton.Value == true
                    tracksFinal = tracksFinal([tracksFinal.pole] == 0);
                end
                if app.zsectioningCheckBox.Value == true
                    tracksFinal = tracksFinal([tracksFinal.z] == 0);
                end
                %calculating D, alpha, MSD
                tracksFinal = calc_D(tracksFinal, frameTime, twenty_five_percent);    %calculating D values
                
                %making D-histogram
                D_list = [tracksFinal.D];     %list of D-values
                figure;
                D_hist = histogram(D_list);
                title('Cytoplasm D-histogram');
                xlabel('D (µm^2/s)');
                ylabel('counts');
                savefig('cytoplasm_D_hist.fig');       %saving figure of histogram
                exportgraphics(gca,'cytoplasm_D_hist.jpg','Resolution',300);
                cytoplasm_D_hist_load = imread('cytoplasm_D_hist.jpg');
                %edges = D_hist.BinEdges
                %counts = D_hist.BinCounts
                %D_list
                %save('cytoplasm_D_histogram','edges','counts','D_list');     %saving histogram info
                app.Image2.ImageSource = cytoplasm_D_hist_load;        %displaying histogram on GUI
                close
                delete('cytoplasm_D_hist.jpg')
                
                %making alpha-histogram
                Alpha_list = [tracksFinal.Alpha];     %list of alpha-values
                figure;
                Alpha_hist = histogram(Alpha_list);
                title('Cytoplasm α-histogram');
                xlabel('α');
                ylabel('counts');
                savefig('cytoplasm_Alpha_hist.fig');       %saving figure of histogram
                exportgraphics(gca,'cytoplasm_Alpha_hist.jpg','Resolution',300);
                cytoplasm_Alpha_hist_load = imread('cytoplasm_Alpha_hist.jpg');
                %edges = Alpha_hist.BinEdges
                %counts = Alpha_hist.BinCounts
                %save('cytoplasm_Alpha_histogram','edges','counts','D_list');     %saving histogram info
                app.Image5.ImageSource = cytoplasm_Alpha_hist_load;        %displaying histogram on GUI
                close
                delete('cytoplasm_Alpha_hist.jpg')
                
                %making MSD graph
                MSD_list = [tracksFinal.MSD];    %list of MSD for each step (averaged of all tracks)
                time = linspace(frameTime,frameTime.*length(MSD_list(:,1)),length(MSD_list(:,1)));
                save('MSD_list','MSD_list')
                figure;
                if app.linearButton.Value == true
                    plot(time,MSD_list);
                end
                if app.loglogButton.Value==true
                    loglog(time,MSD_list);
                end
                title('Mean Squared Displacement plot');
                xlabel('Time (s)');
                ylabel('MSD (µm^2)')
                savefig('cytoplasm_MSD_plot.fig')
                exportgraphics(gca,'cytoplasm_MSD_plot.jpg','Resolution',300);
                cytoplasm_MSD_plot_load = imread('cytoplasm_MSD_plot.jpg');
                app.Image4.ImageSource = cytoplasm_MSD_plot_load;
                close
                delete('cytoplasm_MSD_plot.jpg')
                
                %making probability plot of xNorm from first value of each track
                save('xNorm_values','xNorm_values');
                [mids, h] = myHist2(xNorm_values, binW);
                figure;
                plot(mids,h);
                title('Cytoplasm xNorm Plot')
                xlabel('x-position (µm)')
                ylabel('probability')
                savefig('cytoplasm_xNorm.fig');
                exportgraphics(gca,'cytoplasm_xNorm.jpg','Resolution',300);
                close
                %displaying xNorm plot
                cytoplasm_xNorm_load = imread('cytoplasm_xNorm.jpg');
                app.Image.ImageSource = cytoplasm_xNorm_load;
                delete('cytoplasm_xNorm.jpg')
                
                save('D_list','D_list')
                save('Alpha_list','Alpha_list')
                save('tracksFinal','tracksFinal')
                save('data','tracksFinal','D_list','Alpha_list','MSD_list','xNorm_values') %save the variables for future use
                save('parameters','parameters')
                delete('*_list.mat')
                delete('tracksFinal.mat')
                delete('xNorm_values.mat')
            end
            
            %Membrane Simulation
            if app.MembraneButton.Value == true
                type = 'Membrane';
                %defining variables from user inputs
                %The length unit is micrometer [um] and time unit is seconds [s]                
                R = app.CellRadiusmEditField.Value;   %cell radius
                L = app.CellLengthmEditField.Value;   %cell length
                N = app.offramesEditField.Value;  %number of frames
                frameTime = app.frametimemsEditField.Value/1000;   %frame time
                nTracks = app.oftracksEditField.Value;         %number of tracks
                loc_err = app.localizationerrormEditField.Value;   %localization error
                D = app.Dm2sEditField.Value;  %diffusion coefficient
                a = app.EditField.Value;  %anomalous diffusion coefficient
                binW = app.BinWidthEditField.Value;    %xNorm graph bin width
                percent = "100";      %percent of particles on the membrane
                motion_blur_timing = floor(frameTime.*1000/3)/1000;    %motion blurring timing
                T_exp = N.*frameTime;      %total time
                twenty_five_percent = floor(N/4);          %25 percent of the number of frames (for D and Alpha calculations)
                if app.MotionBlurringCheckBox.Value == true
                    avg = 1;
                end
                if app.MotionBlurringCheckBox.Value == false          %disabling motion blurring if user chooses
                    avg = 0;
                end
                z_min = app.z_minEditField.Value;
                z_max = app.z_maxEditField.Value;
                parameters = struct('Type',type,'frame_time',frameTime,'motion_blur_microstep', motion_blur_timing,'nFrames',N,'total_time',T_exp,'nTracks',nTracks,'L',L,'R',R,'motion_blur',avg,'D',D,'loc_err',loc_err,'z_min',z_min,'z_max',z_max);
                
                %getting xNorm plot
                %xNorm_plots_runner('R',R,'L', L, 'N', nTracks,'factor',1.5,'percent',percent,'locErr',loc_err,'binW',binW);
                %xNorm_load = imread('xNorm.png');
                %app.Image.ImageSource = xNorm_load;
                %close
                
                %running membrane simulations (with nTracks number of tracks) and putting the output D-values into a histogram
                n_tracks = zeros(1,nTracks);
                for i = 1:nTracks
                    n_tracks(i) = nTracks;
                end
                
                clear tracksFinal
                
                tracksFinal = struct([]);
                xNorm_values = []
                cc_param = struct('R',R,'L',L-2.*R);
                plot_model(cc_param);
                figure(1)
                
                for ii = 1:length(n_tracks)
                    [track,params] = surface_random_walk('D',D,'T',T_exp,'avg',avg,'R',R,'sigma', loc_err, 'L', L-2.*R, 'dt_out', frameTime, 'dt', motion_blur_timing, 'exposureTime', frameTime); %sigma is our localization error (um)
                    track(:,[2 3]) = track(:,[3 2]);                   
                    tracksFinal(ii).pole = any(abs(track(:,2))>=L/2-R);
                    tracksFinal(ii).z = any(track(:,3)>=z_max) + any(track(:,3)<=z_min);
                    if app.excludepolesButton.Value == true
                        if app.zsectioningCheckBox.Value == true
                            if tracksFinal(ii).pole == 0
                                if tracksFinal(ii).z == 0
                                    tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                                    tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                                    xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                                    figure(1)
                                    %hold on
                                    plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                                    daspect([1 1 1])
                                    set(gcf,'color','w');
                                    set(gca,'linewidth',2)
                                    hold on
                                end
                            end
                        end
                        if app.zsectioningCheckBox.Value == false
                            if tracksFinal(ii).pole == 0
                                tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                                tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                                xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                                figure(1)
                                %hold on
                                plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                                daspect([1 1 1])
                                set(gcf,'color','w');
                                set(gca,'linewidth',2)
                                hold on
                            end
                        end
                    end
                    if app.wholecellButton.Value == true
                        if app.zsectioningCheckBox.Value == true
                            if tracksFinal(ii).z == 0
                                tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                                tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                                xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                                figure(1)
                                %hold on
                                plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                                daspect([1 1 1])
                                set(gcf,'color','w');
                                set(gca,'linewidth',2)
                                hold on
                            end
                        end
                        if app.zsectioningCheckBox.Value == false
                            tracksFinal(ii).tracksCoordXYZ = track; %with localization error
                            tracksFinal(ii).tracksCoordXY = track(:,1:2:end); %This is because of the nameing, we want to have the projected along cell minor axis and along cell major axis
                            xNorm_values = [xNorm_values tracksFinal(ii).tracksCoordXY(:,1)];      %getting the x-values of each track for the xNorm plot
                            figure(1)
                            %hold on
                            plot3(track(:,1),track(:,2),track(:,3),'.-','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize', 5)     %plotting simulated tracks
                            daspect([1 1 1])
                            set(gcf,'color','w');
                            set(gca,'linewidth',2)
                            hold on
                        end
                    end
                end
                camup([0 1 0])
                xlabel('x (µm)')
                ylabel('y (µm)')
                zlabel('z (µm)')
                savefig('membrane_tracks.fig');
                exportgraphics(gca,'membrane_tracks.jpg','Resolution',300);
                close
                %displaying simulated tracks
                membrane_tracks_load = imread('membrane_tracks.jpg');
                app.Image3.ImageSource = membrane_tracks_load;
                delete('membrane_tracks.jpg')
                
                if app.excludepolesButton.Value == true
                    tracksFinal = tracksFinal([tracksFinal.pole] == 0);
                end
                if app.zsectioningCheckBox.Value == true
                    tracksFinal = tracksFinal([tracksFinal.z] == 0);
                end
                %calculating D, alpha, MSD
                tracksFinal=calc_D(tracksFinal, frameTime, twenty_five_percent); %this is basically calculating the D, alpha, etc from the track points, this also includes a comparison for the 3D values
                %entire cell 2D          
                
                %making D-histogram
                D_list = [tracksFinal.D];    %list of D-values
                figure;
                D_hist = histogram(D_list);
                title('Membrane D-histogram');
                xlabel('D (µm^2/s)');
                ylabel('counts');
                savefig('membrane_D_hist.fig');       %saving figure of histogram
                exportgraphics(gca,'membrane_D_hist.jpg','Resolution',300);
                membrane_D_hist_load = imread('membrane_D_hist.jpg');
                %edges = D_hist.BinEdges
                %counts = D_hist.BinCounts
                %D_list
                %save('membrane_D_histogram','edges','counts','D_list');     %saving histogram info
                app.Image2.ImageSource = membrane_D_hist_load;        %displaying histogram on GUI
                close
                delete('membrane_D_hist.jpg')
                
                %making alpha-histogram
                Alpha_list = [tracksFinal.Alpha];     %list of D-values
                figure;
                Alpha_hist = histogram(Alpha_list);
                title('Membrane α-histogram');
                xlabel('α');
                ylabel('counts');
                savefig('membrane_Alpha_hist.fig');       %saving figure of histogram
                exportgraphics(gca,'membrane_Alpha_hist.jpg','Resolution',300);
                membrane_Alpha_hist_load = imread('membrane_Alpha_hist.jpg');
                %edges = Alpha_hist.BinEdges
                %counts = Alpha_hist.BinCounts
                %save('membrane_Alpha_histogram','edges','counts','D_list');     %saving histogram info
                app.Image5.ImageSource = membrane_Alpha_hist_load;        %displaying histogram on GUI
                close
                delete('membrane_Alpha_hist.jpg')
                
                %making MSD graph
                MSD_list = [tracksFinal.MSD];    %list of MSD for each step (averaged of all tracks)
                time = linspace(frameTime,frameTime.*length(MSD_list(:,1)),length(MSD_list(:,1)));
                save('MSD_list','MSD_list')
                figure;
                if app.linearButton.Value == true
                    plot(time,MSD_list);
                end
                if app.loglogButton.Value==true
                    loglog(time,MSD_list);
                end
                title('Mean Squared Displacement plot');
                xlabel('Time (s)');
                ylabel('MSD (µm^2)')
                savefig('membrane_MSD_plot.fig')
                exportgraphics(gca,'membrane_MSD_plot.jpg','Resolution',300);
                membrane_MSD_plot_load = imread('membrane_MSD_plot.jpg');
                app.Image4.ImageSource = membrane_MSD_plot_load;
                close
                delete('membrane_MSD_plot.jpg')
                
                %making probability plot of xNorm from first value of each track
                save('xNorm_values','xNorm_values')
                [mids, h] = myHist2(xNorm_values, binW);
                figure;
                plot(mids,h);
                title('Membrane xNorm Plot')
                xlabel('x-position (µm)')
                ylabel('probability')
                savefig('membrane_xNorm.fig');
                exportgraphics(gca,'membrane_xNorm.jpg','Resolution',300);
                close
                %displaying xNorm plot
                membrane_xNorm_load = imread('membrane_xNorm.jpg');
                app.Image.ImageSource = membrane_xNorm_load;
                delete('membrane_xNorm.jpg')
                
                save('D_list','D_list')
                save('Alpha_list','Alpha_list')
                save('tracksFinal','tracksFinal')
                save('data','tracksFinal','D_list','Alpha_list','MSD_list','xNorm_values') %save the variables for future use
                save('parameters','parameters')
                delete('*_list.mat')
                delete('tracksFinal.mat')
                delete('xNorm_values.mat')
            end
            currDate = strrep(datestr(datetime), ':', '_');
            mkdir(currDate);
            movefile('data.mat',string(currDate));
            movefile('parameters.mat',string(currDate));
            movefile('*.fig',string(currDate));
            save('currDate','currDate')
        end

        % Button pushed function: UpdateButton_2
        function UpdateButton_2Pushed(app, event)
            binW = app.BinWidthEditField.Value;
            load('currDate')
            addpath(string(currDate))
            load('data.mat')
            %making probability plot of xNorm from first value of each track
            [mids, h] = myHist2(xNorm_values, binW);
            figure;
            plot(mids,h);
            title('xNorm')
            xlabel('x-position')
            ylabel('probability')
            if app.MembraneButton.Value == true
                savefig('membrane_xNorm.fig');
                exportgraphics(gca,'membrane_xNorm.jpg','Resolution',300);
                close
                %displaying xNorm plot
                membrane_xNorm_load = imread('membrane_xNorm.jpg');
                app.Image.ImageSource = membrane_xNorm_load;
            end
            if app.CytoplasmButton.Value == true
                savefig('cytoplasm_xNorm.fig');
                exportgraphics(gca,'cytoplasm_xNorm.jpg','Resolution',300);
                close
                %displaying xNorm plot
                cytoplasm_xNorm_load = imread('cytoplasm_xNorm.jpg');
                app.Image.ImageSource = cytoplasm_xNorm_load;
            end
            delete('*.jpg')
            movefile('*.fig',string(currDate));
            clear
        end

        % Button pushed function: UpdateButton
        function UpdateButtonPushed(app, event)
            frameTime = app.frametimemsEditField.Value;
            load('currDate')
            addpath(string(currDate))
            load('data.mat')
            if app.loglogButton.Value == true
                time = linspace(frameTime,frameTime.*length(MSD_list(:,1)),length(MSD_list(:,1)));
                figure;
                loglog(time,MSD_list);
                title('Mean Squared Displacement plot');
                xlabel('Time (s)');
                ylabel('MSD (µm^2)')
            end
            if app.linearButton.Value == true
                time = linspace(frameTime,frameTime.*length(MSD_list(:,1)),length(MSD_list(:,1)));
                figure;
                plot(time,MSD_list);
                title('Mean Squared Displacement plot');
                xlabel('Time (s)');
                ylabel('MSD')
            end
            if app.MembraneButton.Value == true
                savefig('membrane_MSD_plot.fig')
                exportgraphics(gca,'membrane_MSD_plot.jpg','Resolution',300);
                membrane_MSD_plot_load = imread('membrane_MSD_plot.jpg');
                app.Image4.ImageSource = membrane_MSD_plot_load;
                delete('membrane_MSD_plot.jpg')
            end
            if app.CytoplasmButton.Value == true
                savefig('cytoplasm_MSD_plot.fig')
                exportgraphics(gca,'cytoplasm_MSD_plot.jpg','Resolution',300);
                cytoplasm_MSD_plot_load = imread('cytoplasm_MSD_plot.jpg');
                app.Image4.ImageSource = cytoplasm_MSD_plot_load;
                delete('cytoplasm_MSD_plot.jpg')
            end
            movefile('*.fig',string(currDate));
            close
            clear
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 823 617];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [322 15 487 600];

            % Create simulatedcellwithtracksTab
            app.simulatedcellwithtracksTab = uitab(app.TabGroup);
            app.simulatedcellwithtracksTab.Title = 'simulated cell with tracks';
            app.simulatedcellwithtracksTab.BackgroundColor = [1 1 1];

            % Create Image3
            app.Image3 = uiimage(app.simulatedcellwithtracksTab);
            app.Image3.Position = [20 40 446 520];

            % Create MSDvstTab
            app.MSDvstTab = uitab(app.TabGroup);
            app.MSDvstTab.Title = 'MSD vs. t';
            app.MSDvstTab.BackgroundColor = [1 1 1];

            % Create PlottingOptionsButtonGroup
            app.PlottingOptionsButtonGroup = uibuttongroup(app.MSDvstTab);
            app.PlottingOptionsButtonGroup.Title = 'Plotting Options';
            app.PlottingOptionsButtonGroup.Position = [185 83 123 73];

            % Create loglogButton
            app.loglogButton = uiradiobutton(app.PlottingOptionsButtonGroup);
            app.loglogButton.Text = 'log-log';
            app.loglogButton.Position = [11 27 58 22];
            app.loglogButton.Value = true;

            % Create linearButton
            app.linearButton = uiradiobutton(app.PlottingOptionsButtonGroup);
            app.linearButton.Text = 'linear';
            app.linearButton.Position = [11 5 65 22];

            % Create UpdateButton
            app.UpdateButton = uibutton(app.MSDvstTab, 'push');
            app.UpdateButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateButtonPushed, true);
            app.UpdateButton.BackgroundColor = [0.8 0.8 0.8];
            app.UpdateButton.Position = [197 62 100 22];
            app.UpdateButton.Text = 'Update';

            % Create Image4
            app.Image4 = uiimage(app.MSDvstTab);
            app.Image4.Position = [-85 170 663 390];

            % Create DhistogramsTab
            app.DhistogramsTab = uitab(app.TabGroup);
            app.DhistogramsTab.Title = 'D, α histograms';
            app.DhistogramsTab.BackgroundColor = [1 1 1];

            % Create Image2
            app.Image2 = uiimage(app.DhistogramsTab);
            app.Image2.Position = [66 291 361 279];

            % Create Image5
            app.Image5 = uiimage(app.DhistogramsTab);
            app.Image5.Position = [44 1 406 280];

            % Create xNormTab
            app.xNormTab = uitab(app.TabGroup);
            app.xNormTab.Title = 'xNorm';
            app.xNormTab.BackgroundColor = [1 1 1];

            % Create Image
            app.Image = uiimage(app.xNormTab);
            app.Image.Position = [38 201 437 359];

            % Create BinWidthEditField
            app.BinWidthEditField = uieditfield(app.xNormTab, 'numeric');
            app.BinWidthEditField.Position = [233 162 100 22];
            app.BinWidthEditField.Value = 0.05;

            % Create BinWidthEditFieldLabel
            app.BinWidthEditFieldLabel = uilabel(app.xNormTab);
            app.BinWidthEditFieldLabel.HorizontalAlignment = 'right';
            app.BinWidthEditFieldLabel.Position = [161 162 57 22];
            app.BinWidthEditFieldLabel.Text = 'Bin Width';

            % Create UpdateButton_2
            app.UpdateButton_2 = uibutton(app.xNormTab, 'push');
            app.UpdateButton_2.ButtonPushedFcn = createCallbackFcn(app, @UpdateButton_2Pushed, true);
            app.UpdateButton_2.Position = [193 134 100 22];
            app.UpdateButton_2.Text = 'Update';

            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.UIFigure);
            app.ButtonGroup.AutoResizeChildren = 'off';
            app.ButtonGroup.Title = 'Simulation Type';
            app.ButtonGroup.BackgroundColor = [1 1 1];
            app.ButtonGroup.Position = [1 539 108 79];

            % Create CytoplasmButton
            app.CytoplasmButton = uiradiobutton(app.ButtonGroup);
            app.CytoplasmButton.Text = 'Cytoplasm';
            app.CytoplasmButton.Position = [11 27 103 22];
            app.CytoplasmButton.Value = true;

            % Create MembraneButton
            app.MembraneButton = uiradiobutton(app.ButtonGroup);
            app.MembraneButton.Text = 'Membrane';
            app.MembraneButton.Position = [11 5 103 22];

            % Create RunButton
            app.RunButton = uibutton(app.UIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.BackgroundColor = [0.8 0.8 0.8];
            app.RunButton.Position = [126 593 100 22];
            app.RunButton.Text = 'Run';

            % Create CellRadiusmEditFieldLabel
            app.CellRadiusmEditFieldLabel = uilabel(app.UIFigure);
            app.CellRadiusmEditFieldLabel.HorizontalAlignment = 'right';
            app.CellRadiusmEditFieldLabel.Position = [14 456 95 22];
            app.CellRadiusmEditFieldLabel.Text = 'Cell Radius (µm)';

            % Create CellRadiusmEditField
            app.CellRadiusmEditField = uieditfield(app.UIFigure, 'numeric');
            app.CellRadiusmEditField.Position = [124 456 100 22];
            app.CellRadiusmEditField.Value = 0.5;

            % Create CellLengthmEditFieldLabel
            app.CellLengthmEditFieldLabel = uilabel(app.UIFigure);
            app.CellLengthmEditFieldLabel.HorizontalAlignment = 'right';
            app.CellLengthmEditFieldLabel.Position = [13 419 94 22];
            app.CellLengthmEditFieldLabel.Text = 'Cell Length (µm)';

            % Create CellLengthmEditField
            app.CellLengthmEditField = uieditfield(app.UIFigure, 'numeric');
            app.CellLengthmEditField.Position = [122 419 100 22];
            app.CellLengthmEditField.Value = 3;

            % Create MotionBlurringCheckBox
            app.MotionBlurringCheckBox = uicheckbox(app.UIFigure);
            app.MotionBlurringCheckBox.Text = 'Motion Blurring';
            app.MotionBlurringCheckBox.Position = [22 493 103 22];

            % Create offramesEditFieldLabel
            app.offramesEditFieldLabel = uilabel(app.UIFigure);
            app.offramesEditFieldLabel.HorizontalAlignment = 'right';
            app.offramesEditFieldLabel.Position = [12 385 66 22];
            app.offramesEditFieldLabel.Text = '# of frames';

            % Create offramesEditField
            app.offramesEditField = uieditfield(app.UIFigure, 'numeric');
            app.offramesEditField.Position = [93 385 100 22];
            app.offramesEditField.Value = 20;

            % Create frametimemsEditFieldLabel
            app.frametimemsEditFieldLabel = uilabel(app.UIFigure);
            app.frametimemsEditFieldLabel.HorizontalAlignment = 'right';
            app.frametimemsEditFieldLabel.Position = [10 352 89 22];
            app.frametimemsEditFieldLabel.Text = 'frame time (ms)';

            % Create frametimemsEditField
            app.frametimemsEditField = uieditfield(app.UIFigure, 'numeric');
            app.frametimemsEditField.Position = [114 352 100 22];
            app.frametimemsEditField.Value = 21.7;

            % Create Dm2sEditFieldLabel
            app.Dm2sEditFieldLabel = uilabel(app.UIFigure);
            app.Dm2sEditFieldLabel.HorizontalAlignment = 'right';
            app.Dm2sEditFieldLabel.Position = [13 245 64 22];
            app.Dm2sEditFieldLabel.Text = 'D (µm^2/s)';

            % Create Dm2sEditField
            app.Dm2sEditField = uieditfield(app.UIFigure, 'numeric');
            app.Dm2sEditField.Position = [92 245 100 22];
            app.Dm2sEditField.Value = 0.5;

            % Create EditFieldLabel
            app.EditFieldLabel = uilabel(app.UIFigure);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Position = [3 210 25 22];
            app.EditFieldLabel.Text = 'α ';

            % Create EditField
            app.EditField = uieditfield(app.UIFigure, 'numeric');
            app.EditField.Editable = 'off';
            app.EditField.Position = [43 210 100 22];
            app.EditField.Value = 1;

            % Create localizationerrormEditFieldLabel
            app.localizationerrormEditFieldLabel = uilabel(app.UIFigure);
            app.localizationerrormEditFieldLabel.HorizontalAlignment = 'right';
            app.localizationerrormEditFieldLabel.Position = [12 280 122 22];
            app.localizationerrormEditFieldLabel.Text = 'localization error (µm)';

            % Create localizationerrormEditField
            app.localizationerrormEditField = uieditfield(app.UIFigure, 'numeric');
            app.localizationerrormEditField.Position = [149 280 100 22];

            % Create ButtonGroup_2
            app.ButtonGroup_2 = uibuttongroup(app.UIFigure);
            app.ButtonGroup_2.AutoResizeChildren = 'off';
            app.ButtonGroup_2.Title = 'Analysis Options';
            app.ButtonGroup_2.BackgroundColor = [1 1 1];
            app.ButtonGroup_2.Position = [183 509 114 79];

            % Create wholecellButton
            app.wholecellButton = uiradiobutton(app.ButtonGroup_2);
            app.wholecellButton.Text = 'whole cell';
            app.wholecellButton.Position = [11 27 103 22];
            app.wholecellButton.Value = true;

            % Create excludepolesButton
            app.excludepolesButton = uiradiobutton(app.ButtonGroup_2);
            app.excludepolesButton.Text = 'exclude poles';
            app.excludepolesButton.Position = [11 5 103 22];

            % Create oftracksEditFieldLabel
            app.oftracksEditFieldLabel = uilabel(app.UIFigure);
            app.oftracksEditFieldLabel.HorizontalAlignment = 'right';
            app.oftracksEditFieldLabel.Position = [13 316 61 22];
            app.oftracksEditFieldLabel.Text = '# of tracks';

            % Create oftracksEditField
            app.oftracksEditField = uieditfield(app.UIFigure, 'numeric');
            app.oftracksEditField.Position = [89 316 100 22];
            app.oftracksEditField.Value = 50;

            % Create z_minEditFieldLabel
            app.z_minEditFieldLabel = uilabel(app.UIFigure);
            app.z_minEditFieldLabel.HorizontalAlignment = 'right';
            app.z_minEditFieldLabel.Position = [9 149 38 22];
            app.z_minEditFieldLabel.Text = 'z_min';

            % Create z_minEditField
            app.z_minEditField = uieditfield(app.UIFigure, 'numeric');
            app.z_minEditField.Position = [62 149 100 22];
            app.z_minEditField.Value = -1;

            % Create z_maxEditFieldLabel
            app.z_maxEditFieldLabel = uilabel(app.UIFigure);
            app.z_maxEditFieldLabel.HorizontalAlignment = 'right';
            app.z_maxEditFieldLabel.Position = [9 122 41 22];
            app.z_maxEditFieldLabel.Text = 'z_max';

            % Create z_maxEditField
            app.z_maxEditField = uieditfield(app.UIFigure, 'numeric');
            app.z_maxEditField.Position = [65 122 100 22];
            app.z_maxEditField.Value = 1;

            % Create zsectioningCheckBox
            app.zsectioningCheckBox = uicheckbox(app.UIFigure);
            app.zsectioningCheckBox.Text = 'z-sectioning';
            app.zsectioningCheckBox.Position = [14 177 86 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Diffusion_Simulator_Original
         

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end