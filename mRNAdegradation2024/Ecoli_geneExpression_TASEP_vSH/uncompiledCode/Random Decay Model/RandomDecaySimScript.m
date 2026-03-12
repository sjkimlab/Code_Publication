%--------------------------------------------------------------------------------------------------
% Code written by: Albur Hassan
% sjkimlab at University of Illinois at Urbana Champaign.
% Original creation date: 2026/2/26
% Last edited: 2026/03/09
% Contact: ahassan4@illinois.edu
%--------------------------------------------------------------------------------------------------
% Program Description: script to run the simulations for the Random Decay situation.
% This script runs the Random Decay simulation. in order to change a parameter in the simulation 
% add in a line of code following the param = struct; line such as:
%      param.kRiboLoading= 0.025;
% The default paramaters for the script can be found in the Ecoli_geneExpression_TASEP_function.m 
% file. This script should save a file with the name {date}_{time}__geneExpressionTASEP.mat containing
% the details of the simulation
%--------------------------------------------------------------------------------------------------


%clear off the final_output variable sincce that will be used to store data.
clear final_output;
%param.ribofollowing = 0 is needed to set the random decay model.
%ideally you can set the param file to be the final_output.parameters variable from a previous run and change the paramaters needed.
param = struct;
param.ribofollowing = 0;
nloci= 10; %%% 10 will take about 10 min.
final_output= Ecoli_geneExpression_TASEP_function(param);
fish1_signals = zeros(size(final_output.fishTime,1),nloci,2);
fish1_signals(:,1,1) = final_output.fish5Signal;
fish1_signals(:,1,2) = final_output.fish3Signal;
final_output.allProteinSignals=zeros(size(final_output.fishTime,1),nloci);
final_output.allProteinSignals(:,1)= final_output.proteinSignal;

SEndSEQOutput = struct;
SEndSEQOutput(1).RNAP_exitTimes = final_output.RNAP_exitTimes;
SEndSEQOutput(1).deg_exitTimes = final_output.deg_exitTimes;
for i=2:nloci
	output = Ecoli_geneExpression_TASEP_function(param);
    fish1_signals(:,i,1) = output.fish5Signal;
    fish1_signals(:,i,2) = output.fish3Signal;
    final_output.allProteinSignals(:,i)= output.proteinSignal;
    final_output.proteinSignal = final_output.proteinSignal + output.proteinSignal;
    final_output.transcribing_5signal = final_output.transcribing_5signal +output.transcribing_5signal;
    final_output.DNAStats(:,i) = output.DNAStats;
    final_output.FirstRiboRNAPdistance = [final_output.FirstRiboRNAPdistance;output.FirstRiboRNAPdistance];
    final_output.RiboperRNAP = [final_output.RiboperRNAP;output.RiboperRNAP];
    final_output.RNAP_startTimes = [final_output.RNAP_startTimes,output.RNAP_startTimes];
    final_output.RNAP_endTimes = [final_output.RNAP_endTimes,output.RNAP_endTimes];
    final_output.RNAlifetimes = [final_output.RNAlifetimes,output.RNAlifetimes];
    final_output.prematureterm = final_output.prematureterm + output.prematureterm;
    final_output.RNAlifetimes_analysis =[final_output.RNAlifetimes_analysis,output.RNAlifetimes_analysis];
    final_output.mu1lifetimes =[final_output.mu1lifetimes,output.mu1lifetimes];
    final_output.mu2lifetimes =[final_output.mu2lifetimes,output.mu2lifetimes];
    final_output.ribo_loadTimes = [final_output.ribo_loadTimes;output.ribo_loadTimes];
    final_output.ribo_loadTimes2 = [final_output.ribo_loadTimes2;output.ribo_loadTimes2];
    final_output.RNAseq = final_output.RNAseq + output.RNAseq;
    final_output.crazyfish = final_output.crazyfish +output.crazyfish;
    final_output.NETseq = final_output.NETseq + output.NETseq;
    final_output.Sendseq = [final_output.Sendseq;output.Sendseq];
    if mod(i, 10) == 0
    	SEndSEQOutput(floor((i+1)/10)).RNAP_exitTimes = output.RNAP_exitTimes;
    	SEndSEQOutput(floor((i+1)/10)).deg_exitTimes = output.deg_exitTimes;
    end
end
final_output.crazyfish = final_output.crazyfish/nloci;
final_output.fish5Signal = mean (fish1_signals(:,:,1),2);
final_output.fish5FanoFactor = var(fish1_signals(:,:,1),0,2)./final_output.fish5Signal;
final_output.fish3Signal = mean (fish1_signals(:,:,2),2);
final_output.fish3FanoFactor = var(fish1_signals(:,:,2),0,2)./final_output.fish3Signal;
final_output.fish5Signalerr=std(fish1_signals(:,:,1),0,2)/sqrt(nloci);
final_output.fish3Signalerr=std(fish1_signals(:,:,2),0,2)/sqrt(nloci);
final_output.proteinSignal = final_output.proteinSignal/nloci;
final_output.RNAseq = final_output.RNAseq/nloci;
final_output.NETseq = final_output.NETseq/nloci;
final_output.nloci = nloci;
final_output.SEndSEQOutput =SEndSEQOutput;
filename ="" + string(datetime('now','Format','yyyyMMdd_HHmm'))+"_geneExpressionTASEP";
save((filename+".mat"),'final_output');