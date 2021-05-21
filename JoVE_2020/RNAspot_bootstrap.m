function Xans = RNAspot_bootstrap(dayList,daily_scale,daily_thresh,daily_affine,timecourseInfo)
%{
-About-
This code calculates mRNA numbers per cell based on total spot intensities. Spot magnitude is normalized
by a spot intensity of a single mRNA and spots within a cell is summed.
This code do bootstrap by sampling an equal number of cells from different
day's experiments.

Note: parfor can cause an issue depending on the machine. See comments with
%%%%%%% mark.

-Inputs-
dayList: contains file path of each day
daily_scale: single mRNA normalization factor, different for each day.
daily_thresh: threshold for spot curation, different for each day.
daily_affine: affine transformation matrix for each day
timecouseInfo: actual time (minutes) of FISH sample; day and the file name
for the images from that time.

-Outputs-
Xans: Bootstrap result of spot properties. 
Result is also saved as mat file.

-Example-
Xans = RNApercell_by_spotmag_bootstrap(dayList,daily_scale,daily_thresh,timecourseInfo)
see FISHworkflow for the input
  
-Supplementary-

-Keywords-
spotmag, SM FISH analysis, bootstrap

-Dependencies-
FISHworkflow

-References-

-Author-
Sangjin Kim, modified 2021, May 10
%}

timePoints = timecourseInfo.time;
days = timecourseInfo.day;
files = timecourseInfo.file;
SMnormalization = 1; % 1= do normalization in spotThresholding function. Choose 0 to avoid normalization

% Go over each time point
for ti = 1:length(timePoints)
    cellindex = []; 
    cellNum = [];

    positionArray1a = []; dArray1a = []; lArray1a = []; magnitudeArray1a = []; 
    stepAreaArray1a = []; stepLengthArray1a = []; xArray1a = []; yArray1a = [];
    affine1a1 = []; affine2a1 = [];
    celllengthArray1 = [];
    
    positionArray1b = []; dArray1b = []; lArray1b = []; magnitudeArray1b = []; 
    stepAreaArray1b = []; stepLengthArray1b = []; xArray1b = []; yArray1b = [];
    affine1b1 = []; affine2b1 = [];
    
    
    for di = 1:length(days{ti})
        dd = days{ti}(di);
        dayfolder = dayList(dd).folderpath;
        scaleZ5 = daily_scale(dd,1); scaleZ3 = daily_scale(dd,2);
        thresh1 = daily_thresh(dd,1); thresh2 = daily_thresh(dd,2);
        ABaffine1 = daily_affine(dd).ABaffine1;
        ABaffine2 = daily_affine(dd).ABaffine2;
        ACaffine1 = daily_affine(dd).ACaffine1;
        ACaffine2 = daily_affine(dd).ACaffine2;
        
        tifile = strcat('t',sprintf('%01.0f',files{ti}(di)),'-mesha');
        cellListfile = strcat(dayfolder,tifile,'Cy5.mat'); load(cellListfile); TEMPcellLista = cellList;
        cellListfile = strcat(dayfolder,tifile,'Cy3.mat'); load(cellListfile); TEMPcellListb = cellList;
        TEMPcellindex = CellCounter(TEMPcellLista,dd); %take [day, frame, cellNum] of non-empty cell in the cellList
        cellindex{dd} = TEMPcellindex(:,2:3);
        cellNum(dd) = length(TEMPcellindex); % total 'real' cell # in the cellList

        if SMnormalization == 0
            scaleZ5 = 1; scaleZ3 = 1;
        end;
        
        %add affine matrix to cellList. For example, add ABaffine (YFP-Cy5)
        %to Cy5 cellList.
        TEMPcellLista = addAffinetoCellList(TEMPcellLista,ABaffine1,ABaffine2); 
        TEMPcellListb = addAffinetoCellList(TEMPcellListb,ACaffine1,ACaffine2);
        
        [TEMPcelllengthArraya,TEMPpositionArraya,TEMPdArraya,TEMPlArraya,TEMPmagnitudeArraya,...
         TEMPstepAreaArraya,TEMPstepLengthArraya,TEMPxArraya,TEMPyArraya,TEMPaffine1a,TEMPaffine2a] = spotThresholding(TEMPcellLista,1,thresh1,scaleZ5);
        [TEMPcelllengthArrayb,TEMPpositionArrayb,TEMPdArrayb,TEMPlArrayb,TEMPmagnitudeArrayb,...
         TEMPstepAreaArrayb,TEMPstepLengthArrayb,TEMPxArrayb,TEMPyArrayb,TEMPaffine1b,TEMPaffine2b] = spotThresholding(TEMPcellListb,2,thresh2,scaleZ3);

        %make arrays of spot properties. dd = day (for a given t) a= Cy5 channel, b = Cy3 channel.
        celllengthArray1{dd} = TEMPcelllengthArraya;
        positionArray1a{dd} = TEMPpositionArraya; positionArray1b{dd} = TEMPpositionArrayb;
        dArray1a{dd} = TEMPdArraya; dArray1b{dd} = TEMPdArrayb;
        lArray1a{dd} = TEMPlArraya; lArray1b{dd} = TEMPlArrayb;
        magnitudeArray1a{dd} = TEMPmagnitudeArraya; magnitudeArray1b{dd} = TEMPmagnitudeArrayb;
        stepAreaArray1a{dd} = TEMPstepAreaArraya; stepAreaArray1b{dd} = TEMPstepAreaArrayb;
        stepLengthArray1a{dd} = TEMPstepLengthArraya; stepLengthArray1b{dd} = TEMPstepLengthArrayb;
        xArray1a{dd} = TEMPxArraya; xArray1b{dd} = TEMPxArrayb;
        yArray1a{dd} = TEMPyArraya; yArray1b{dd} = TEMPyArrayb;
        affine1a1{dd} = TEMPaffine1a; affine1b1{dd} = TEMPaffine1b;
        affine2a1{dd} = TEMPaffine2a; affine2b1{dd} = TEMPaffine2b;
        
        clear TEMP*;
    end;
   
    close all;
    nboot = min(cellNum(find(cellNum>0))); %equal number of cells withdrawn each day
        
    numberOfIterations = 3000; %for bootstrapping
    
    bootspot1CHa = zeros(numberOfIterations,9); %bootstrap result for Cy5. Change the size according to the output of spotStatCH1.
    bootspot1CHb = zeros(numberOfIterations,9); %bootstrap result for Cy3. Change the size according to the output of spotStatCH1. 
    bootspot2CHa = zeros(numberOfIterations,11); %bootstrap result for Cy5-Cy3
    bootspot2CHb = zeros(numberOfIterations,11); %bootstrap result for Cy3-Cy5
    
    parallelobject = parpool(12); %%%%%%%%%% comment this if regular for is used instead of parfor
    parfor i = 1:numberOfIterations   %change parfor to for if debugging
        %makeTestArrays randomly samples nboot number of cells from each
        %day and generates a temperary array.
        %The temporary array is analyzed in spotStat1CH and spotStat2CH.
        %parfor is to repeat this random sampling and analysis
        [celllengthArray, positionArraya,dArraya,lArraya,magnitudeArraya,...
        stepAreaArraya,stepLengthArraya,xArraya,yArraya,...
        affine1a,affine2a,positionArrayb,dArrayb,lArrayb,magnitudeArrayb,...
        stepAreaArrayb,stepLengthArrayb,xArrayb,yArrayb,...
        affine1b,affine2b] = makeTestArrays(nboot,cellindex,celllengthArray1, positionArray1a,dArray1a,lArray1a,magnitudeArray1a,stepAreaArray1a,...
                                            xArray1a,stepLengthArray1a,yArray1a,affine1a1,affine2a1,...
                                            positionArray1b,dArray1b,lArray1b,magnitudeArray1b,stepAreaArray1b,...
                                            xArray1b,stepLengthArray1b,yArray1b,affine1b1,affine2b1);
        %output from makeTestArrays is a list of {cellNum} from different
        %days and frames.
        bootspot1CHa(i,:) = spotStat1CH(celllengthArray{1}, positionArraya{1},dArraya{1},lArraya{1},...
                                        magnitudeArraya{1},stepAreaArraya{1},stepLengthArraya{1});
        bootspot1CHb(i,:) = spotStat1CH(celllengthArray{1}, positionArrayb{1},dArrayb{1},lArrayb{1},...
                                        magnitudeArrayb{1},stepAreaArrayb{1},stepLengthArrayb{1});
        bootspot2CHa(i,:) = spotStat2CH(xArraya,yArraya,magnitudeArraya,xArrayb,...
                                        yArrayb,magnitudeArrayb,affine1a,affine2a,affine1b,affine2b);
        bootspot2CHb(i,:) = spotStat2CH(xArrayb,yArrayb,magnitudeArrayb,xArraya,...
                                        yArraya,magnitudeArraya,affine1b,affine2b,affine1a,affine2a);
    end
    delete(parallelobject); %%%%%%%%%% comment this if regular for is used instead of parfor
    
    %cell-level statistics-------------------------------
    %CH1 sum of spot mag per cell = RNA per cell
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,7))),7); 
    cell_ans(ti,1) = mean(temp); cell_ans(ti,2) = std(temp);
    
    %CH2 sum of spot mag per cell = RNA per cell
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,7))),7); 
    cell_ans(ti,3) = mean(temp); cell_ans(ti,4) = std(temp);
    
    %total number of cells analyzed (nboot=equal number of cells withdrawn each day)
    cell_ans(ti,5) = nboot*length(find(cellNum>0)); 
    
    %spot-level statistics-------------------------
    %Channel 1 (CH1) only (Cy5)
    %CH1 spot intensity (mean and ste)
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,1))),1); 
    spot1_ans(ti,1) = mean(temp); spot1_ans(ti,2) = std(temp);
    
    %CH1 spot intensity (std)
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,2))),2); 
    spot1_ans(ti,3) = mean(temp); 
    
    %CH1 cells with spot (%)
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,4))),4);
    spot1_ans(ti,4) = mean(temp); spot1_ans(ti,5) = std(temp);
    
    %CH1 spot # per cell
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,5))),5); 
    spot1_ans(ti,6) = mean(temp); spot1_ans(ti,7) = std(temp);
    
    %CH1 spot # per cell (>0, among cells with at least one spot)
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,6))),6); 
    spot1_ans(ti,8) = mean(temp); spot1_ans(ti,9) = std(temp);
    
    %CH1 avg localizationD (mean and ste)
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,8))),8); 
    spot1_ans(ti,10) = mean(temp); spot1_ans(ti,11) = std(temp);
    
    %CH1 avg localizationD (std)
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,9))),9);
    spot1_ans(ti,12) = mean(temp);
    
    %CH1 total # spot analyzed per bootstrap 
    temp = bootspot1CHa(find(~isnan(bootspot1CHa(:,3))),3); 
    spot1_ans(ti,13) = mean(temp); spot1_ans(ti,14) = std(temp);
    
%     %CH1 spot localization D
%     %c = (0.05:0.1:1.05); %localization along short axis (d)
%     for j = 15:25, 
%         spot1_ans(ti,j) = mean(bootspot1CHa(find(~isnan(bootspot1CHa(:,(j-5)))),(j-5)));
%     end;
       
    %spot-level statistics-------------------------
    %Channel 2 only (Cy3)
    %CH2 spot intensity (mean and ste)
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,1))),1); 
    spot2_ans(ti,1) = mean(temp); spot2_ans(ti,2) = std(temp);
    
    %CH2 spot intensity (std)
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,2))),2); 
    spot2_ans(ti,3) = mean(temp); 
    
    %CH2 cells with spot (%)
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,4))),4);
    spot2_ans(ti,4) = mean(temp); spot2_ans(ti,5) = std(temp);
    
    %CH2 spot # per cell
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,5))),5); 
    spot2_ans(ti,6) = mean(temp); spot2_ans(ti,7) = std(temp);
    
    %CH2 spot # per cell (>0, among cells with at least one spot)
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,6))),6); 
    spot2_ans(ti,8) = mean(temp); spot2_ans(ti,9) = std(temp);
    
    %CH2 avg localizationD (mean and ste)
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,8))),8); 
    spot2_ans(ti,10) = mean(temp); spot2_ans(ti,11) = std(temp);
    
    %CH2 avg localizationD (std)
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,9))),9);
    spot2_ans(ti,12) = mean(temp);
    
    %CH2 total # spot analyzed per bootstrap 
    temp = bootspot1CHb(find(~isnan(bootspot1CHb(:,3))),3); 
    spot2_ans(ti,13) = mean(temp); spot2_ans(ti,14) = std(temp);
    
%     %CH2 spot localization D
%     %c = (0.05:0.1:1.05); %localization along short axis (d)
%     for j = 15:25, 
%         spot2_ans(ti,j) = mean(bootspot1CHb(find(~isnan(bootspot1CHb(:,(j-5)))),(j-5)));
%     end;    
    
    %%%%%%%%%%%%%%%%%
    %spot-level statistics-------------------------
    %Co-localization analysis: For a given spot in CH1 (or A channel), is there a CH2 (or B channel) spot co-localized?
    %Percentage of spotA in a cell without any B spot
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,1))),1); % (%) of type 1
    spot12_ans(ti,1) = mean(temp); spot12_ans(ti,2) = std(temp);
    
    %Percentage of spotA uncolocalized with any Bspot (even if there is Bspot in the cell)
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,2))),2); % (%) of type 2
    spot12_ans(ti,3) = mean(temp); spot12_ans(ti,4) = std(temp);
    
    %Percentage of spotA colocalized with Bspot
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,3))),3); % (%) of type 3(coloc ratio)
    spot12_ans(ti,5) = mean(temp); spot12_ans(ti,6) = std(temp);
    
    %Intensity of spotA colocalized with Bspot
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,4))),4); % coloc spot int
    spot12_ans(ti,7) = mean(temp); spot12_ans(ti,8) = std(temp);

    %Ratio of spot intensity A/B
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,6))),6); % coloc spot ratio Ch1/ch2
    spot12_ans(ti,9) = mean(temp); spot12_ans(ti,10) = std(temp);
    
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,8))),8); % coloc spot num per cell1
    spot12_ans(ti,11) = mean(temp); spot12_ans(ti,12) = std(temp);
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,9))),9); % uncoloc spot num per cell1
    spot12_ans(ti,13) = mean(temp); spot12_ans(ti,14) = std(temp);
    
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,10))),10); % coloc spot num per cell among all cells
    spot12_ans(ti,15) = mean(temp); spot12_ans(ti,16) = std(temp);
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,11))),11); % uncoloc spot num per cell  among all cells
    spot12_ans(ti,17) = mean(temp); spot12_ans(ti,18) = std(temp);
    
    %%%%%%%%%%%%%%%%%
    %spot-level statistics-------------------------
    %Co-localization analysis: For a given spot in CH2 (or B channel), is there a CH1 (or A channel) spot co-localized?
    %Percentage of spotB in a cell without any A spot
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,1))),1); % (%) of type 1
    spot21_ans(ti,1) = mean(temp); spot21_ans(ti,2) = std(temp);
    
    %Percentage of spotB uncolocalized with any Aspot (even if there is Aspot in the cell)
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,2))),2); % (%) of type 2
    spot21_ans(ti,3) = mean(temp); spot21_ans(ti,4) = std(temp);
    
    %Percentage of spotB colocalized with Aspot
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,3))),3); % (%) of type 3(coloc ratio)
    spot21_ans(ti,5) = mean(temp); spot21_ans(ti,6) = std(temp);
    
    %Intensity of spotB colocalized with Aspot
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,4))),4); % coloc spot int
    spot21_ans(ti,7) = mean(temp); spot21_ans(ti,8) = std(temp);
    
    %Ratio of spot intensity B/A
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,6))),6); % coloc spot ratio Ch2/ch1
    spot21_ans(ti,9) = mean(temp); spot21_ans(ti,10) = std(temp);
 
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,8))),8); % coloc spot num per cell1
    spot21_ans(ti,11) = mean(temp); spot21_ans(ti,12) = std(temp);
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,9))),9); % uncoloc spot num per cell1
    spot21_ans(ti,13) = mean(temp); spot21_ans(ti,14) = std(temp);
    
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,10))),10); % coloc spot num per cell among all cells
    spot21_ans(ti,15) = mean(temp); spot21_ans(ti,16) = std(temp);
    temp = bootspot2CHa(find(~isnan(bootspot2CHa(:,11))),11); % uncoloc spot num per cell  among all cells
    spot21_ans(ti,17) = mean(temp); spot21_ans(ti,18) = std(temp);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    display(ti)
end

timePoints= timePoints';
Xans = [timePoints, cell_ans, spot1_ans, spot2_ans, spot12_ans, spot21_ans];

save(['RNAspot_bootstrap_result.mat'],'timePoints','dayList','cell_ans','spot1_ans','spot2_ans','spot12_ans','spot21_ans');

end % for this function


function cellindex = CellCounter(cellList,n)
i = 0;
for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>4
            i = i+1;
            cellindex(i,1) = n;
            cellindex(i,2) = frame;
            cellindex(i,3) = cell;
        end
    end
end
end %for CellCounter

function cellList = addAffinetoCellList(cellList,A,B)
for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>4
            cellList{frame}{cell}.affine1 = A;
            cellList{frame}{cell}.affine2 = B;
        end
    end
end
end %for addAffinetoCellList

function [celllengthArray,positionArray,dArray,lArray,magnitudeArray,...
          stepAreaArray,stepLengthArray,xArray,yArray,affine1,affine2] = spotThresholding(cellList,ch,thresh,scale)
%output arrays are in {frame}{cellNum} format.

sizeOfCellList  = cellfun(@length,cellList);
celllengthArray = cell(1,length(cellList)); 
positionArray   = cell(1,length(cellList)); 
dArray          = cell(1,length(cellList)); 
lArray          = cell(1,length(cellList)); 
magnitudeArray  = cell(1,length(cellList)); 
stepAreaArray   = cell(1,length(cellList)); 
stepLengthArray = cell(1,length(cellList)); 
xArray          = cell(1,length(cellList)); 
yArray          = cell(1,length(cellList)); 
affine1         = cell(1,length(cellList));
affine2         = cell(1,length(cellList));
for frame=1:length(cellList)
    celllengthArray{frame} = cell(1,sizeOfCellList(frame));
    positionArray{frame}   = cell(1,sizeOfCellList(frame));
    dArray{frame}          = cell(1,sizeOfCellList(frame));
    lArray{frame}          = cell(1,sizeOfCellList(frame));
    magnitudeArray{frame}  = cell(1,sizeOfCellList(frame));
    stepAreaArray{frame}   = cell(1,sizeOfCellList(frame));
    stepLengthArray{frame} = cell(1,sizeOfCellList(frame));
    xArray{frame}          = cell(1,sizeOfCellList(frame));
    yArray{frame}          = cell(1,sizeOfCellList(frame));
    affine1{frame}         = cell(1,sizeOfCellList(frame));
    affine2{frame}         = cell(1,sizeOfCellList(frame));
    for cellNum=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cellNum}) && length(cellList{frame}{cellNum}.mesh)>4 && isfield(cellList{frame}{cellNum},'spots') 
            
            %---for given cell, spot threshold
            if ch == 1
                prf = cellList{frame}{cellNum}.signal1./cellList{frame}{cellNum}.steparea;
            elseif ch == 2
                prf = cellList{frame}{cellNum}.signal2./cellList{frame}{cellNum}.steparea;
            end;
            
            for i=1:3
                prf = prf*0.5+prf([1 1:end-1])*0.25+prf([2:end end])*0.25; %smoothing
            end
            
            pos = cellList{frame}{cellNum}.spots.positions;
            prf = [0;prf];
            
            ind = prf(pos+1)>thresh; %thresholding (remove weak spots): 1 or 0
            ind2 = cellList{frame}{cellNum}.spots.magnitude>0 & cellList{frame}{cellNum}.spots.magnitude<20;
            ind2 = ind2';
            try
                ind = ind & ind2; %logic numbers
            catch
            end;

            celllengthArray{frame}{cellNum} = cellList{frame}{cellNum}.length;
            positionArray{frame}{cellNum} = cellList{frame}{cellNum}.spots.positions(ind); 
            xArray{frame}{cellNum} = cellList{frame}{cellNum}.spots.x(ind);
            yArray{frame}{cellNum} = cellList{frame}{cellNum}.spots.y(ind);
            lArray{frame}{cellNum} = cellList{frame}{cellNum}.spots.l(ind)/cellList{frame}{cellNum}.length;
            %spot position along the long axis. normalized by the cell length!
            dArray{frame}{cellNum} = cellList{frame}{cellNum}.spots.d(ind);
            magnitudeArray{frame}{cellNum} = cellList{frame}{cellNum}.spots.magnitude(ind)/scale; 
            %spot intensity, normalization by SM, or scale
            stepAreaArray{frame}{cellNum}   = cellList{frame}{cellNum}.steparea;
            stepLengthArray{frame}{cellNum} = cellList{frame}{cellNum}.steplength;
            affine1{frame}{cellNum}         = cellList{frame}{cellNum}.affine1;
            affine2{frame}{cellNum}         = cellList{frame}{cellNum}.affine2;
        end;
    end;
end;
end %for spotThresholding


function [celllengthArray,positionArraya,dArraya,lArraya,magnitudeArraya,...
    stepAreaArraya,stepLengthArraya,xArraya,yArraya,...
    affine1a,affine2a,positionArrayb,dArrayb,lArrayb,magnitudeArrayb,...
    stepAreaArrayb,stepLengthArrayb,xArrayb,yArrayb,...
    affine1b,affine2b] = makeTestArrays(nboot,cellindex,celllengthArray1,positionArray1a,dArray1a,lArray1a,magnitudeArray1a,stepAreaArray1a,...
                                            xArray1a,stepLengthArray1a,yArray1a,affine1a1,affine2a1,...
                                            positionArray1b,dArray1b,lArray1b,magnitudeArray1b,stepAreaArray1b,...
                                            xArray1b,stepLengthArray1b,yArray1b,affine1b1,affine2b1)
%input arrays are in {day}{frame}{cellNum} format.
%output arrays are in {cellNum} format. 
%i.e., day, frame information are lost in the output. Output is a series of cells from different
%days and frames.

for di = 1:length(cellindex)
    try
        cellNumperday(di) = length(cellindex{di});
    catch
    end;
end;
% nboot = min(cellNumperday(find(cellNumperday>0)));

testcellindex = [];
for di = 1:length(cellNumperday) %go over each day
    if cellNumperday(di)>0
        %randomly choose nboot cells from a given day
        sampleindx = unidrnd(cellNumperday(di),nboot,1); %random sampling
        sampleindx2 = zeros(nboot,3);
        sampleindx2(:,1) = di; % for day
        sampleindx2(:,2:3) = cellindex{di}(sampleindx,:); %for frame and cell
        testcellindex = [testcellindex; sampleindx2];
    end;
end;
    
for j = 1:size(testcellindex,1) %generate temporary list of spot properties
    celllengthArray{1}{j}   = celllengthArray1{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    positionArraya{1}{j}   = positionArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    dArraya{1}{j}          = dArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    lArraya{1}{j}          = lArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    magnitudeArraya{1}{j}  = magnitudeArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    stepAreaArraya{1}{j}   = stepAreaArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    stepLengthArraya{1}{j} = stepLengthArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    xArraya{1}{j}          = xArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    yArraya{1}{j}          = yArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    affine1a{1}{j}         = affine1a1{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    affine2a{1}{j}         = affine2a1{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    positionArrayb{1}{j}   = positionArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    dArrayb{1}{j}          = dArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    lArrayb{1}{j}          = lArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    magnitudeArrayb{1}{j}  = magnitudeArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    stepAreaArrayb{1}{j}   = stepAreaArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    stepLengthArrayb{1}{j} = stepLengthArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    xArrayb{1}{j}          = xArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    yArrayb{1}{j}          = yArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    affine1b{1}{j}         = affine1b1{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    affine2b{1}{j}         = affine2b1{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
end;
end %for makeTestArrays

function Xans = spotStat1CH(celllengthArray, positionArray,dArray,lArray,magnitudeArray,stepAreaArray,stepLengthArray)
%The input arrays are a series of cells (not spots).
%They are converted to a series of spots in this function (i.e. spotList).

spotArrayNonEmpty = find(~cellfun('isempty',magnitudeArray)>0); %indexes of magnitudeArray
tempMagnitudeArray = magnitudeArray(spotArrayNonEmpty);
tempMagnitudeArray = cell2mat(tempMagnitudeArray);

if length(tempMagnitudeArray)>0
    
    spotList = zeros(length(tempMagnitudeArray),3); %length of spots, or the number of spots
    spotList(:,1) = tempMagnitudeArray; %spotList(:,1) spot magnitude

    spotList(:,3) = cell2mat(lArray(spotArrayNonEmpty));  %spotList(:,3) normalized L
    %spot position in long axis (L)
    %already normalized by cell length in spotThresholding

    positionArrayZeroValueIndex = cell2mat(cellfun(@(x) x==0,positionArray(spotArrayNonEmpty),'uniformOutput',0));
    try  
        spotList(positionArrayZeroValueIndex,2)= 1;
    catch
    end;
    
    %spotList(:,2) normalized D. Normalization is done by the segment width where the spot is located. (stepArea/stepLength)
    dArrayInPositionArray      = cell2mat(dArray(spotArrayNonEmpty)); 
    stepAreaInPositionArray    = stepAreaArray(spotArrayNonEmpty); %still cell array
    stepLengthInPositionArray  = stepLengthArray(spotArrayNonEmpty); 
    tempPositionArray          = positionArray(spotArrayNonEmpty);
    counter = 0;
    for ii = 1:length(stepAreaInPositionArray) %ii = cell
        localPositionArray = tempPositionArray{ii};
        for jj = 1:length(localPositionArray)
            counter = counter + 1;
            try %% incase there is a spot with position (seg num) = 0
            tempValue = stepAreaInPositionArray{ii}(localPositionArray(jj))./stepLengthInPositionArray{ii}(localPositionArray(jj))/2;
            catch 
                continue; %skip the remaining commands in the for loop
            end
            spotList(counter,2)= abs(dArrayInPositionArray(counter)/tempValue); %spotList(:,2) normalied d
        end
    end;
    
    %sum spot magnitudes in each cell
    totalmagpercell = cell2mat(cellfun(@sum,magnitudeArray,'uniformoutput',0));
    %count number of spots per cell
    spotspercell = cell2mat(cellfun(@length,positionArray,'uniformoutput',0));
    %divide the number of cells with spot(s) by the total number of cells
    cellWithSpot = length(spotArrayNonEmpty)*100/length(positionArray); %%%%%%%%

    %calculate some numbers as an output. 
    %one can change this section to output different quantities and to be
    %added to the bootstrap statistics in the main function.
    Xans(1) = mean(spotList(:,1)); % mean spot intensity
    Xans(2) = std(spotList(:,1)); % std of spot intensity
    Xans(3) = length(spotList); % total spot number
    Xans(4) = cellWithSpot;  %cells with spot %
    Xans(5) = mean(spotspercell); % spot # per cell including zeros
    Xans(6) = mean(spotspercell(spotspercell>0)); % spot # per cell excluding zeros
    Xans(7) = mean(totalmagpercell); %sum of spot magnitude per cell = RNA per cell
    Xans(8) = mean(spotList(:,2)); % mean of spot localization (d)
    Xans(9) = std(spotList(:,2)); %std of spot localization (d)
    %---example: histograms can be calculated.
%     c = (0.05:0.1:1.05); %localization
%     Xans(10:20) = hist(spotList(:,2),c)*100/length(spotList(:,1)); %localization along short axis (d)
%     %figure, bar(c,Xans(7:17));
%     binc =(0:1:5);
%     Xans(21:26) = hist(spotList(:,1),binc)*100/length(spotList(:,1));

else % if there is no spot
    Xans(1) = NaN; Xans(2) = NaN; Xans(3) = 0; Xans(4) = 0; Xans(5) = 0; Xans(6) = NaN;
    Xans(7) = 0; Xans(8) = NaN; Xans(9) = NaN;    
%     Xans(10:20) = NaN;     
%     Xans(21:26) = NaN;
end;
end %for spotStat1CH


function Xans = spotStat2CH(xArraya,yArraya,magnitudeArraya,xArrayb,yArrayb,magnitudeArrayb,...
                            affine1a,affine2a,affine1b,affine2b)
%The input arrays are a series of cells (not spots).
%Spots in A channel is analyzed for the presence and absence of
%co-localization (with a spot in B channel) --> see "AspotList".

colocThresh =150/64.2; %criteria for spot colocalization (pixels) 150 nm= camera pixel size.

spotAix = 0;
cellcount = 0; cellAB = 0; cellA = 0; AspotList = [];  
for cellNum=1:size(xArraya{1},2)
    cellcount = cellcount + 1; 
    colocstate = 0; %flag for colocstate for a given spot (of A channel)
    colocnumCell2(cellcount) = 0; %colocalized spot number (A channel) per cell (among all cells)
    uncolocnumCell2(cellcount) = 0; %colocalized spot number (A channel) per cell (among all cells)
    
    if ~isempty(xArraya{1}{cellNum})
        cellA = cellA + 1; % cell number with a spot in A channel
        Aspot.x0 = xArraya{1}{cellNum}; Aspot.y0 = yArraya{1}{cellNum}; Aspot.mag = magnitudeArraya{1}{cellNum};
        colocnumCell1(cellA) = 0; %set colocalized spot number = 0; unless found later.
        uncolocnumCell1(cellA) = 0; %set uncolocalized spot num = 0;
        
        %check if there is a spot in B channel.
        if ~isempty(xArrayb{1}{cellNum})
            cellAB = cellAB + 1; %cell number with spots in both channels
%                     colocnumCell2(AB) = 0; % colocalized spot number = 0; unless found later.
            %save the original (x,y) coorindate of spots in A channel.
            Aspot.xy0 = [Aspot.x0',Aspot.y0',ones(size(Aspot.x0,2),1)];
            %affine transform (x,y) coordinate of spots in A channel.
            Aspot.x = Aspot.xy0*affine1a{1}{cellNum}; Aspot.x = Aspot.x';
            Aspot.y = Aspot.xy0*affine2a{1}{cellNum}; Aspot.y = Aspot.y';
            %affine transform (x,y) coordinate of spots in B channel.
            Bspot.x0 = xArrayb{1}{cellNum}; Bspot.y0 = yArrayb{1}{cellNum}; Bspot.mag = magnitudeArrayb{1}{cellNum};
            Bspot.xy0 = [Bspot.x0',Bspot.y0',ones(size(Bspot.x0,2),1)];
            Bspot.x = Bspot.xy0*affine1b{1}{cellNum}; Bspot.x = Bspot.x';
            Bspot.y = Bspot.xy0*affine2b{1}{cellNum}; Bspot.y = Bspot.y';

            %for a given A-channel spot (i), check if there is colocalized
            %spot in B channel.
            for i = 1:length(Aspot.x)
                L = []; clear tmp*;
                for j = 1:length(Bspot.x)
                    L = [L; sqrt((Aspot.x(i)-Bspot.x(j))^2+(Aspot.y(i)-Bspot.y(j))^2)];
                end;
                [tmpL,tmpi] = min(L);

                spotAix = spotAix+1; %for each A spot
                if tmpL < colocThresh %criteria for spot colocalization
                    colocstate = 2; % if A and B colocalized
                    AspotList(spotAix,:) = [colocstate,Aspot.mag(i), Bspot.mag(tmpi)]; %save magnitude of colocalized spots in A and B channels.
                    colocnumCell1(cellA) = colocnumCell1(cellA) + 1; %colocalized spot number per cell (among cells with A spot)
                    colocnumCell2(cellcount) = colocnumCell2(cellcount) + 1; %colocalized spot number per cell (among all cells)
                else %if A spot do not colocalize with B
                    colocstate = 1; %if A and B didn't colocalize
                    AspotList(spotAix,:) = [colocstate,Aspot.mag(i),0];
                end;
            end;
        else % if there is no B spot in the cell
            colocstate = 0.5; %if only A exist
            for i = 1:length(Aspot.x0)
                spotAix = spotAix+1; 
                AspotList(spotAix,:) = [colocstate,Aspot.mag(i),0];
            end;

        end;
        uncolocnumCell1(cellA) = length(Aspot.x0) - colocnumCell1(cellA); %number of A spots that do not have colocalized counter part in B channel within a cell.
        uncolocnumCell2(cellcount) = length(Aspot.x0) -colocnumCell2(cellcount); %number of A spots that do not have colocalized counter part in B channel within a cell.
    end;
            
end

if length(AspotList)>0 
    %percent of A spot in a cell without B spot.
    Xans(1) = length(find(AspotList(:,1) == 0.5))*100/length(AspotList);  %colocstate = 0.5; only A exist
    %percent of A spot not colocalized with B spot.
    Xans(2) = length(find(AspotList(:,1) == 1))*100/length(AspotList);  %colocstate = 1; B exist in the cell, but do not coloc with A;
    %percent of A spot colocalized with a B spot.
    Xans(3) = length(find(AspotList(:,1) == 2))*100/length(AspotList);  %colocstate = 2; colocalized
    %sum of Xans(1)+Xans(2)+Xans(3) should be 100.
    
    Xans(4) = mean(AspotList(find(AspotList(:,1) == 2),2)); %mean of colocalized A signal intensity;
    Xans(5) = std(AspotList(find(AspotList(:,1) == 2),2));
    
    %mean(AspotList(find(AspotList(:,1)~= 2),2)); %mean of uncolocalized A signal intensity;
    
    tmp = AspotList(find(AspotList(:,1) == 2),2)./AspotList(find(AspotList(:,1) == 2),3);
    Xans(6) = mean(tmp); %ratio of intensities (A/B)
    Xans(7) = std(tmp);
    %-Histogram version:
    % binc =(0:1:5);
    % Xans(7:12) = hist(tmp,binc)*100/length(tmp);

    Xans(8) = mean(colocnumCell1); % colocalized spot number per cell (among cells with A-channel spot)
    Xans(9) = mean(uncolocnumCell1); % 'alone' spot number per cell (among cells with A-channel spot)
    Xans(10) = mean(colocnumCell2); % colocalized spot number per cell (among all cells)
    Xans(11) = mean(uncolocnumCell2); % 'alone' spot number per cell (among all cells)
    %-Histogram version:
%     binc = (0:1:5);
%     Xans(14:19) = hist(colocnumCell2,binc')*100/length(colocnumCell2);
%     Xans(20:25) = hist(uncolocnumCell2,binc')*100/length(uncolocnumCell2);
else
    Xans(1:11) = NaN;
end;
end %for spot2CH
