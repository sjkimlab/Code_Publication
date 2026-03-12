t = 1500;
all_RNA= [];RNA_STATE=[];
geneLength = final_output.parameters.geneLength;
SEndSEQOutput = final_output.SEndSEQOutput;
for i =1:30 %length(SEndSEQOutput)
	temp = SEndSEQOutput(i).RNAP_exitTimes<t & SEndSEQOutput(i).deg_exitTimes>t;
	temp_RNA_STATE= max(SEndSEQOutput(i).RNAP_exitTimes(:,:))>t;
	temp_RNA_STATE= temp_RNA_STATE + 2*(sum(SEndSEQOutput(i).RNAP_exitTimes(:,:)>0,1)<geneLength &max(SEndSEQOutput(i).RNAP_exitTimes(:,:))<t);
	temp_RNA_STATE= temp_RNA_STATE + 3*(SEndSEQOutput(i).RNAP_exitTimes(end,:)<t& SEndSEQOutput(i).RNAP_exitTimes(end,:)>0);
	condition = sum(temp>0,1)>0;
	all_RNA = [all_RNA,temp(:,condition)];RNA_STATE = [RNA_STATE,temp_RNA_STATE(condition)];
end


figure('position',[0,0,450,600]);
pbaspect([1 1 1]); 
set(gca,'FontSize', 16);
hold on;
xlabel('Gene position (bp)');
ylabel('RNA #');
subplot(4,1,1);
pbaspect([4 1 1]); 
xlim([0,geneLength-1]);
ylim([0,17]);
hold on;
plot(final_output.RNAseq(:,1));
ylabel('RNA #');
subplot(4,1,2:4);
pbaspect([1 1 1]); 
hold on;
xlabel('Gene position (bp)');
ylabel('RNA #');
%yticks([0,2000,4000,6000]);
xlim([0,geneLength-1]);
ylim([0,1025]);
xticks([0:1000:3000]);
yticks([0:200:1000]);
title("sendSEQ data");

%further trim the data
steps= 50;
condition = sum(all_RNA(1:steps:end,:)>0,1)>0;

% Example signal
%all_RNA_grouped = movsum(all_RNA, [0 49], 1);
%all_RNA_grouped =all_RNA_grouped(1:steps:geneLength,:);
%all_RNA_grouped= all_RNA_grouped./steps;
%condition = sum(all_RNA_grouped(:,:)>=0.5,1)>0;

all_RNA = all_RNA(:,condition);RNA_STATE= RNA_STATE(condition);
%lets sort first
[temp,idx] = sort(sum(all_RNA,1));
RNA_STATE= RNA_STATE(idx);all_RNA= all_RNA(:,idx);

[RNA_STATE,idx] = sort(RNA_STATE);
all_RNA= all_RNA(:,idx);

for i = 1:length(RNA_STATE)
	for j = 1:steps:geneLength
		if(all_RNA(j,i)>= 0.5)  %_grouped(ceil(j/steps),i)
			x = [j , j+steps ,j+steps ,j];
			y= [i-1, i-1, i ,i];
			if(RNA_STATE(i) ==1)
				C = [0.2 0.7 0.2];
			elseif(RNA_STATE(i) ==2)
				C= [0.9 0.45 0.1];
			else 
				C = [0.9 0.1 0.1];
			end
			fill(x,y,C,'EdgeColor','none','FaceAlpha',0.7);
		end
	end
end

%xlim([0,max(sendSEQ(:,2))]);
%ylim([0,NascentRNA+ReleasedRNA]);

%{
totalRNA = size(sendSEQ,1);
for i = 1:totalRNA
	x = [sendSEQ(i,1) , sendSEQ(i,2) ,sendSEQ(i,2) ,sendSEQ(i,1)];
	y= [i-1, i-1, i ,i];
	fill(x,y,'r','EdgeColor','none','FaceAlpha',0.7);
end
%}
