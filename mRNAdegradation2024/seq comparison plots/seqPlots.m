%{
---------------------------------------------------------------------------
Author: Sangjin Kim - sangjin@illinois.edu
    Creation date: 10/20/2024 
    Last updated at 12/17/2025

Description: this script is to plot SEnd-seq PF data and Rif-seq halflife
data for the paper.

====== Requirements =======
Make sure Source Data 3 and 4 are downloaded from the paper and placed in the same
folder as this script.

====== Output =======
Figures

10/20/2024: shixinPlot1
11/03/2024: shixinBCMdata4
11/30/2024: YZanalysis1, YZgene3
12/17/2025: combine and clean up
---------------------------------------------------------------------------
%}

%% 1. Plot TE and PF correlation
%% TE and PF n = 1336 processed by Xiangwu Ju
dataB1=readtable('Source Data 3.xlsx','Sheet',2,'Range','B:C','ReadVariableNames',true); %TE, PF
TE = table2array(dataB1(:,1)); %TE
PF = table2array(dataB1(:,2));%PF

DotDensityMedian(TE,PF,'TE','PF',0.05)


%% 2. Examine BCM effect on PF
%% Control (-BCM) and BCM (+BCM) processed by Xiangwu Ju
dataB2=readtable('Source Data 3.xlsx','Sheet',3,'Range','C:D','ReadVariableNames',true); %PF control, PF+BCM
PF_con = table2array(dataB2(:,1)); %PF without BCM
PF_BCM = table2array(dataB2(:,2)); %PF with BCM

DotDensityMedian(PF_con,PF_BCM,'PF_{-BCM}','PF_{+BCM}',0.05); hold on;
glineX = [0,3];
plot(glineX, glineX,'y-'); hold off;


%Histogram of ratio +BCM/-BCM
lowPF = find(PF_con<1);
highPF = find(PF_con>=1);
BCMdiff = PF_BCM./PF_con; 
figure, 
histogram(BCMdiff(highPF,1),'Normalization','probability','BinEdges',-3:0.1:3 ,'FaceColor','g','FaceAlpha',0.5); hold on;
histogram(BCMdiff(lowPF,1),'Normalization','probability','BinEdges',-3:0.1:3 , 'FaceColor','r','FaceAlpha',0.5);
xlabel('Ratio of PF_{+BCM} over PF_{-BCM}'); ylabel('Prob');
xlim([0 3]); legend({'High PF','Low PF'});
[h,p] = kstest2(BCMdiff(highPF,1),BCMdiff(lowPF,1))


%% 3. Examine mRNA half-life (Rif-seq) and TE (Ribo-seq)
%% data processed by Yan Zhang
%n = 820 one-phase genes (shown in Figure 4b)
dataB3=readtable('Source Data 4.xlsx','Sheet',3,'Range','B:C','ReadVariableNames',true); %halflife, TE
HL_1ph = table2array(dataB3(:,1)); %Half-life of one-phase decay genes
TE_1ph = table2array(dataB3(:,2)); %TE of one-phase decay genes

DotDensityMedian(TE_1ph,HL_1ph,'TE','mRNA half-life (min)',0.05); 
title(sprintf('One-phase decay genes (n = %d)', length(HL_1ph)));
corr(TE_1ph,HL_1ph)
corr(TE_1ph,HL_1ph, 'Type','Spearman')

%n = 123 two-phase genes (shown in Extended Data Figure
dataB4=readtable('Source Data 4.xlsx','Sheet',4,'Range','B:D','ReadVariableNames',true); %halflife, TE
HL_2phA = table2array(dataB4(:,1)); %Early Half-life of two-phase decay genes (short)
HL_2phB = table2array(dataB4(:,2)); %Late Half-life of two-phase decay genes (long)
TE_2ph = table2array(dataB4(:,3)); %TE of one-phase decay genes

%Combine one-phase genes + late half-life of two-phase decay genes
%(Extended Data Figure 7f)
DotDensityMedian([TE_1ph;TE_2ph],[HL_1ph;HL_2phB;],'TE','mRNA half-life (min)',0.05); 
title(sprintf('All genes (n = %d) of late phase', length(HL_1ph)+length(HL_2phB)));
corr([TE_1ph;TE_2ph],[HL_1ph;HL_2phB;])
corr([TE_1ph;TE_2ph],[HL_1ph;HL_2phB;],'Type','Spearman')

%Combine one-phase genes + early half-life of two-phase decay genes
%(Extended Data Figure 7g)
DotDensityMedian([TE_1ph;TE_2ph],[HL_1ph;HL_2phA;],'TE','mRNA half-life (min)',0.05); 
title(sprintf('All genes (n = %d) of early phase', length(HL_1ph)+length(HL_2phA)));
corr([TE_1ph;TE_2ph],[HL_1ph;HL_2phA;])
corr([TE_1ph;TE_2ph],[HL_1ph;HL_2phA;],'Type','Spearman')

%mRNA decay pattern with TE and PF
%Fisher's exact test (Figure 4d) n = 682
dataB5=readtable('Source Data 4.xlsx','Sheet',5,'Range','B:E','ReadVariableNames',true); %HLa, HLb, TE, PF
HLa = table2array(dataB5(:,1)); %Early Half-life 
HLb = table2array(dataB5(:,2)); %Late Half-life
TE = table2array(dataB5(:,3)); %TE 
PF = table2array(dataB5(:,4)); %TE 

medTE = median(TE);
medPF = median(PF);
One_ph = find(HLa == HLb); %index for one-phase genes
Two_ph = find(HLa ~= HLb); %index for two-phase genes

%Create array for fisher test1. Median of TE is 1.11
%Is there nonrandom association between one/two phase decay trait and
%high/low TE?
X1 = length(find(TE(Two_ph)<medTE));
X2 = length(find(TE(Two_ph)>=medTE));
X3 = length(find(TE(One_ph)<medTE));
X4 = length(find(TE(One_ph)>=medTE));
[h,p,stats] = fishertest([X1, X2; X3, X4;])

%Create array for fisher test1. Median of PF is 1.05
%Is there nonrandom association between one/two phase decay trait and
%high/low PF?
X1 = length(find(PF(Two_ph)<medPF));
X2 = length(find(PF(Two_ph)>=medPF));
X3 = length(find(PF(One_ph)<medPF));
X4 = length(find(PF(One_ph)>=medPF));
[h,p,stats] = fishertest([X1, X2; X3, X4;])

