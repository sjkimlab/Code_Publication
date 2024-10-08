%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 9/1/2021
    Last update date: 10/26/2023

10/24/2023: version for combining multiple days of SPT experiments

Description: combineTF combines the tracksFinal files from multiple days of
SPT experiments. It will then calcuate quantities for analysis (e.g. MSD,
Diff, alpha, jump distances). Eventually, these quantities are converted
into matrices (arrays) which are convenient for future plotting


=====Functions=====
    spotNorm_yh.m  (major)
    reNumCells_yh.m
    getCellInfo.m
    getTracksInfo.m
    diffAnalysis.m
    dataToVectors.m
    stepToVectors.m

%}

%% Step 0. Set up the parameters & flags

clear, clc, close all

% flag for combining tracksFinal
%     true: automatically combine tracks in all good (isolated) cells
%     false: manually select cell and only combine tracks in those cells
autoCombineMovieFlag = true; 

% set up the storage path for analysis result, users should set up their own path
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';


%% Step 1. Combine the tracksFinal & cellMesh of multiple movies into one file

tic
clearvars tfComb
cellMeshAll = [];   cellRecord = {};  extraName = ''; 

strainToCombine = string( input( 'Which strains do you want to combine (like: 187 292): ', 's'));
strainNum = split( strainToCombine, ' ')';
strainNum( strainNum=="" ) = []; % in case you input extra space

% find folder of different days
varList = dir( [ varPath 'Variables\SK' char( strainNum) '\2*']);

combNum = 1;
if length( varList) > 1
    fprintf( '\n~~~~~ There are %d days of experiments for this strain ~~~~~\n', length(varList));
    for k = 1:length(varList)
        disp( ['  ' num2str(k) '. ' varList(k).name])
    end
    a = string( input( '\nWhich days do you want to combine(like: 1 3 4): ', 's')); % just input: 1 3 4 5
    combNum = double( split(a, ' ')');
    combNum( isnan(combNum)) = []; % in case you input extra space
end


for j = combNum
    savePath = [varList(j).folder '\' varList(j).name '\'];
            
    list = dir( [savePath '*_Variables.mat']); % list files that match name
    load( [savePath list(1).name]) % load some basic info of the strain
    fprintf( '~~~ Cell & analysis info loaded, start combining tracksFinal ~~~\n\n')
    
    for i = 1: length(list)  % you can also choose a subset of the list
%         savePath = dir( 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\Variables\SK187\');
        
%         savePath = [ varPath 'Variables\' strain '\' expDate extraName '\']; 
        % load _Variables.mat file of individual movie
        clearvars tracksFinal_InIsolatedCell
        load( [savePath list(i).name], 'tracksFinal_InIsolatedCell', 'cellNum',...
            'folderName', 'meshName', 'meshPath', 'cellList', 'cellMesh',...
            'fileName', 'cellFilled', 'goodCells', 'badMovie')    
        fprintf( '  - %s loaded ~~~\n',  fileName);    

        % we only take tracks in isolated cells
        tfGood = tracksFinal_InIsolatedCell;

        % store filename & badMovie flag in new fields of the tfGood structure
        nTracks = length( tfGood);          value = cell( nTracks, 1);    
        value(:) = { string( fileName)};    [tfGood.origin] = value{:};    
        value(:) = { badMovie};             [tfGood.badMovie] = value{:};    

        % store tracksFilled flag (tracks in filled cells) in tfGood
        tracksFilled = num2cell( ismember( cellNum, find( cellFilled)));
        [ tfGood.filled] = tracksFilled{:};

        if ~exist('tfComb', 'var') 
        % create tfComb with same fields as tfGood
            tfComb = tfGood;
            tfComb(:) = [];
        end

        if ~autoCombineMovieFlag % manually select the good cells based on overlay image
            fprintf( '     You can choose from %s\n',  join( string( goodCells), ' '))
            flag = input( '     Are all cells good? (0:Yes, -1:None, [x x]:certain cells): ');
            if flag == -1 % skip this movie
                continue
            elseif isnumeric(flag) % combine tracks from only certain cells
                if ismember( flag, goodCells) % is the input part of goodCells?
                    goodCells = flag;
                    tfGood = tfGood( ismember( cellNum, flag));
                else
                    error('      You have to select from good cells with tracks in it')
                end
            end
        end

        % stack tfGood & cellMesh from different movies together
        tfComb = [tfComb; tfGood]; 
        cellMeshAll = [cellMeshAll; cellMesh( goodCells)]; 
        cellRecord = [cellRecord; {fileName goodCells folderName meshName meshPath}]; % {fileName, goodCells}

        fprintf( '    Included Good cells with tracks: %s\n', join( string( goodCells), ' '))
        if sum( ~cellFilled) > 0
            fprintf( '    Not mostly filled cells: %s\n\n', join( string( find( ~cellFilled')), ' '))    
        else
            fprintf( '    All cells are mostly filled ~~~\n\n')
        end        
    end
end

fprintf( '~~~~~~~ tracksFinal Combination Finished ~~~~~~~\n\n');

%% Step 2. reNumCell & Get cell & tracks info before saving to tracksFinal

Date = input( 'What is the combDate (like ''220728 comb 1 pix''): ', 's');

tracksFinal = reNumCells_yh( tfComb); % reorder the cell number, load Variables

pixelSize = 160e-9;     % Andor camera, unit: m
timeStep = 21.742e-3;   % frame time, unit: s  [check the movie output], needed in MSD fitting

[cellInfo, poleBounds] = getCellInfo( tracksFinal, cellMeshAll, pixelSize);
[nTracks, totalCells, cellNum, tracksLength, tracksOrigin, tracksFilled, badMovie] = getTracksInfo( tracksFinal);

tfPath = [varPath 'Track Analysis Archive\' strain '\'];
if ~exist( tfPath, 'dir')
    mkdir( tfPath)
end

tfName = [strain  ' tracksFinal' extraName ' ' Date];

save( [tfPath tfName], 'tracksFinal', 'cellRecord', 'cellMeshAll',...
    'cellInfo', 'poleBounds', 'nTracks', 'totalCells', 'pixelSize', 'timeStep', ...
    'varPath', 'savePath', 'tfPath', 'tfName', 'Date', 'strain', 'extraName');

toc
fprintf( '~~~ Step 2: ''%s'' saved to **Track Analysis Archive** ~~~\n\n', tfName)


%% Step 3. Diffusion analysis (calculate quantities)

tic

% remove original uTrack output fields in tracksFinal to save space
tracksFinal = rmfield( tracksFinal, {'tracksFeatIndxCG' 'tracksCoordAmpCG' 'seqOfEvents'});

% Calculate MSD & Diff & alpha & Fitting
diffAnalysis

tfAllPath = [ varPath '\tracksFinal\'];

tfAllName = [ 'tracksFinal ' strain extraName ' ' Date ' All']; 
save( [tfAllPath tfAllName], 'tracksFinal', 'cellRecord', 'cellMeshAll',...
    'cellInfo', 'poleBounds', 'nTracks', 'totalCells', 'cellNum', ...
    'tracksLength', 'tracksOrigin', 'tracksFilled', 'badMovie',...
    'pixelSize', 'timeStep', 'fitRange', 'maxTau', 'EnsMSD', 'EnsTAMSD',...
    'varPath', 'savePath', 'tfPath', 'tfName', 'tfAllPath', 'tfAllName',...
    'Date', 'strain', 'extraName');

% toc
fprintf( '~~~ Step 3: ''%s'' saved to varPath ~~~\n\n', tfAllName)


%% Step 4. Convert the tracksFinal structure into Matrix, for plotting

tic
matrixPath = [ varPath '\Matrix\'];
if ~exist( matrixPath, 'dir')
    mkdir( matrixPath)
end

[ tracksLNorm, tracksxNorm, tracksLNorm4, tracksxNorm4, tracksDiff] = dataToVectors( tracksFinal);
% tracksDiff: diff (um^2/s), locErr (nm), alpha, dalpha (um^2/s)

% cell subregion constraints (cell pole or middle part)
bound = poleBounds( cellNum); % cell pole bound for each track
tracksMid = false( nTracks, 2); % based on [first spot, whole track]

LNorm = tracksLNorm; % [first spot, end spot, minLNorm, maxLNorm]

tracksMid(:, 1) = LNorm(:,1) > bound & LNorm(:,1) < 1-bound; % first spot not at pole
tracksMid(:, 2) = LNorm(:,3) > bound & LNorm(:,4) < 1-bound; % whole track not at pole

tracksMid4 = tracksLNorm4 > bound & tracksLNorm4 < 1-bound;

    % get xNorm 12 frame, added 12/26/23
    longT = tracksLength >= 12;
    tf = tracksFinal( longT);

    tracksxNorm12 = nan( length( tf), 12);
    tracksLNorm12 = nan( length( tf), 12);

    for i = 1: length( tf)
        spotNorm = tf(i).spotPosNorm; % [LNorm, xNorm]
        tracksLNorm12(i,:) = spotNorm(1:12, 1);
        tracksxNorm12(i,:) = spotNorm(1:12, 2);
    end

    LNorm = tracksLNorm12;    bound = poleBounds( cellNum( longT));
    inCell = tracksxNorm12 >= -1 & tracksxNorm12 <= 1;
    tracksMid12 = LNorm > bound & LNorm < 1-bound & inCell; % first spot subregion constraint
    
    
matrixName = ['Matrix ' strain extraName ' ' Date];
save( [matrixPath matrixName], 'cellRecord', 'cellMeshAll',...
    'cellInfo', 'poleBounds', 'nTracks', 'totalCells', 'cellNum', ...
    'tracksLength', 'tracksOrigin', 'tracksFilled', 'badMovie',...
    'pixelSize', 'timeStep', 'fitRange', 'maxTau', 'EnsMSD', 'EnsTAMSD',...
    'tracksLNorm', 'tracksxNorm', 'tracksLNorm4', 'tracksxNorm4', 'tracksDiff',...
    'tracksLNorm12', 'tracksxNorm12', 'tracksMid12',... % added 10/26/2023
    'tracksMid', 'tracksMid4', 'folderName', 'Date', 'strain', 'extraName',... 
    'varPath', 'savePath', 'tfPath', 'tfName', 'tfAllPath', 'tfAllName', 'matrixPath', 'matrixName');

toc
fprintf( '~~~ Step 4: ''%s'' saved to ''varPath\\Matrix\\''~~~\n\n', matrixName);


%% Step 5. Calculate the Step Information

tic
stepPath = [varPath 'Steps\'];
if ~exist( stepPath, 'dir')
    mkdir( stepPath)
end

[tracksSteps, tracksStepsLongi, tracksRg, tracksStepsNorm, tracksTheta] = stepToVectors( tracksFinal, cellInfo, pixelSize);

% stack step information from tracks to steps vectors
steps = cell2mat( tracksSteps); % all single step displacement (unit: um)
stepsLongi = cell2mat( tracksStepsLongi);
stepsNorm = cell2mat( tracksStepsNorm); % [LNorm, xNorm] for all steps
stepsTheta = cell2mat( tracksTheta); % theta value for all steps

stepsBound = repelem( bound, tracksLength-1, 1); % replicate the bound for steps in the same track

stepsMid = stepsNorm(:,1) > stepsBound & stepsNorm(:,1) < 1-stepsBound;
stepName = ['Steps ' strain extraName ' ' Date];
save( [stepPath stepName], 'tracksRg', 'tracksMid', ...
    'steps', 'stepsLongi', 'stepsNorm', 'stepsTheta', 'stepsBound', 'stepsMid',...    
    'cellInfo', 'nTracks', 'totalCells', 'cellRecord', 'tracksFilled', 'badMovie',...
    'pixelSize', 'timeStep', 'varPath', 'savePath', 'tfPath', 'tfName',...
    'tfAllPath', 'tfAllName', 'matrixPath', 'matrixName', 'stepPath', 'stepName',...
    'totalCells', 'cellRecord', 'strain', 'Date', 'extraName');

toc
fprintf( '~~~ Step 5: ''%s'' saved to ''varPath\\Steps\\'' ~~~\n\n', stepName);

% tEnd = toc( tStart);
% fprintf( '~~~~~ All Done, Hello World, Total Time: %.1fs ~~~~~~\n\n', tEnd)

 %% Step 6. Display Analysis Result Summary

cellShifted = length( unique( cellNum( ~badMovie)));
cellFilled = length( unique( cellNum( tracksFilled)));

fprintf( ['%s Analysis Result: \n\n     %d total cells with %d tracks, %d shifted cells with %d tracks,\n'...
    '     There are %d shifted filled cells with %d tracks,\n     xNorm (1st spot) has %d tracks not in cap region, in total %d tracks \n'...
    '     xNorm (first 4 spots) has %d tracks not in cap region, in total %d spots\n\n'],...
    folderName, totalCells, nTracks, cellShifted, sum( ~badMovie),...
    cellFilled, sum( tracksFilled), sum( tracksMid(:,1)), sum( ~isnan( tracksxNorm(:,1))),...
    sum( min( tracksMid4, [], 2)), sum( sum( ~isnan( tracksxNorm4( min( tracksMid4, [], 2),:)))))
    
% disp( cellRecord)


