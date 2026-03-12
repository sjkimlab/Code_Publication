function output= Ecoli_geneExpression_TASEP_function(input_parameters)

%--------------------------------------------------------------------------------------------------
% Code written by: Albur Hassan
% sjkimlab at University of Illinois at Urbana Champaign.
% Original creation date: 2024/11/10
% Last edited: 2026/03/05
% Contact: ahassan4@illinois.edu
%--------------------------------------------------------------------------------------------------
% Program Description: Main code for the TASEP simulations. This code runs the TASEP simulation for
% one loci and can be called multiple times to simulate many loci. This function takes in a struct
% which defines its parameters.Below is an example as to how to use this function to simulate a
% loci with kriboloading rate of 0.025 as opposed to its default:
% 
% 	param = struct;
% 	param.kRiboLoading = 0.025;
% 	output = Ecoli_geneExpression_TASEP_function(param);
%
% This function is usually called in dataviewer.mlapp file
%----------------------------------------------

    tStart = tic;
    %default parameters
	parameters =struct;parameters.avgSpeed1 =19;parameters.avgSpeed2 = 19;parameters.ribospeed= 19;
	parameters.degSpeed=19;
	parameters.kLoading =1/20;parameters.kRiboLoading = 0.1;parameters.KRutLoading =0;
	parameters.simtime=2000;parameters.glutime=1600; 
	parameters.mu1=1400;parameters.mu2=120;parameters.muPT=120;parameters.muProtein=6000;
	parameters.ribofollowing = 2;parameters.RNAP_Profile =0;parameters.Random_Model =0;
	parameters.boolRNAPRiboCoupling = 1;parameters.PT_Model = 1; 
	parameters.geneLength=3075;%DNA length 
	parameters.loadOneRibo = 1; parameters.RNAP_dwellTimeProfile = ones(parameters.geneLength,1);
	if nargin >=1 %number of arguments in
		parameterlabels= fieldnames(input_parameters);
		for i = 1:numel(parameterlabels)
				parameters.(parameterlabels{i})= input_parameters.(parameterlabels{i});
		end
	end
	parameters.rutSpeed=5*parameters.ribospeed;



	avgSpeed1=parameters.avgSpeed1;avgSpeed2=parameters.avgSpeed2;riboSpeed = parameters.ribospeed; degSpeed = parameters.degSpeed;
	percentSpeed = 75; %old stuff for the lolz
	geneLength=parameters.geneLength;   
	RNAP_width=35; %in bp
	dx=1;dt=0.1;
	simtime=parameters.simtime;glutime=parameters.glutime; %glu= glucose, sim = simulation 
	nocollision=0;
	mu1=parameters.mu1;mu2=parameters.mu2;muPT=parameters.muPT;muProtein=parameters.muProtein;
	ribofollowing = parameters.ribofollowing; % 0 for non rib following and 1 for last ribo following 2 for 5'->3' ended decay at ribo speed.
	RNAP_Profile =parameters.RNAP_Profile;Random_Model =parameters.Random_Model;
	boolRNAPRiboCoupling = parameters.boolRNAPRiboCoupling;PT_Model = parameters.PT_Model;
	kLoading =parameters.kLoading;kRiboLoading = parameters.kRiboLoading; KRutLoading=parameters.KRutLoading;
    loadOneRibo = parameters.loadOneRibo;

	% statisticTerms[1] =#rna produced statisticTerms[2] =#proteins produced
	statisticTerms =-1*ones(4,1);

	%RiboStatistics[1] = ribo-rnap distances RiboStatistics[2] = ribo per rnap 
	RiboStatistics = -1*ones(2,2);


	% Input parameters for ribosome initiation and elongation
	riboExitTimes =zeros(geneLength/dx,2,2); % lets make it a 3d array such as [bp,RNAP#,Ribo#] %ceil(glutime*kLoading)
	Ribo_locs=[]; % 2d array [RNAP#, Ribo#]
	RNAP_RiboCoupling =[0];% 1 for its coupled 0 for not
	Ribo_width =30;
	rho_width =30;


	%Distribution of mRNA numbers (FISH)
	probe1 = round(660/3075*geneLength); %SH's RT primer ;%40:40:960; %1 % probing 5'-end mRNA
	probe2 = round(2890/3075*geneLength); %SH's RT primer %2000:40:2940;%geneLength; % probing 3'-end mRNA
	fish1_signal =zeros(simtime/dt,1);
	fish2_signal =zeros(simtime/dt,1);
	transcribing_5signal =zeros(simtime/dt,1);
	prematureterm = zeros(simtime/dt,3);% (:,1) total released rna %(:,2) prematurely terminated rna %(:,3) full length rna  
	decayed1 =zeros(simtime/dt,1);
	decayed2 =zeros(simtime/dt,1);


	%rut sites
	%% we dont need to track the exittimes we just need to keep the rho positions
	rut_sites =[round(500*geneLength/3075)];
	rutSpeed=parameters.rutSpeed;
	minRholoadRNA = 80-rho_width;
	rho_locs =[];
	rut_loadT=[];
	specificDwelltimeRho=dx/rutSpeed* ones(geneLength/dx,1);
	tempRho=0;
	r_loc_time=zeros(2,length(rut_sites),2);
	PTpercent =0;
	if PT_Model==1 ||PT_Model==2
		KRutLoading=0;
		PTpercent =parameters.KRutLoading;
	end

	tempExitTimes =[]; %rename to temptimes or something
	%RNAlifetimes=exprnd(mu1);
	RNAlifetimes=[];
	RNAlifetimes_3prime=[];	 % remove all mentions of this to go back to 3' and 5' decay together
	RNAlifetimes_analysis =[];
	avgDwelltime1 = dx/avgSpeed1; %sec per nucleotide
	avgDwelltime2 = dx/avgSpeed2; %sec per nucleotide
	riboavgDwelltime = dx/riboSpeed;degavgDwelltime = dx/degSpeed;
	loadt= exprnd(1/kLoading); %*rand;
	rnap_locs=[];
	Riboloadt= loadt + exprnd(1/kRiboLoading); % 1d array of size sz(RNAP_locs)
	if Random_Model ==1 
		loadt= normrnd(1/kLoading,sqrt(1/kLoading));
		Riboloadt= loadt+  normrnd(1/kRiboLoading,sqrt(1/kRiboLoading));
	end
	RNAlifetimes(1) =loadt+exprnd(mu1); %calculate cootranscriptional lifetime
	RNAlifetimes_3prime(1) =simtime;
	RNAlifetimes_analysis(1)=RNAlifetimes(1) -loadt;
	specificDwelltime1 = avgDwelltime1 * parameters.RNAP_dwellTimeProfile;
	specificDwelltime2 = avgDwelltime2 * parameters.RNAP_dwellTimeProfile;
	RibospecificDwelltime1=riboavgDwelltime*ones(geneLength/dx,1);
	DegspecificDwelltime1 = degavgDwelltime*ones(geneLength/dx,1);
	RNAP_exitTimes = zeros(geneLength/dx,1);
	deg_exitTimes = zeros(geneLength/dx,1);
	deg_locs = [];


	%%time to prep that protein
	Protein_lifetime= [];
	protein_signal=zeros(simtime/dt,1);
	ribo_loadTimes = -ones(2,2);



	%Simulation Loop
	for t = 0:dt:simtime
		%%RNAP loading
		
		if(loadt <=t & t <glutime & ~isempty(rnap_locs) & rnap_locs(length(rnap_locs))-RNAP_width <=0 )  %% make sure no rnap on the site at loading time.
			loadt = t + exprnd(1/kLoading);% calculate load time for next RNAP
		end
		
		if(loadt <=t & t <glutime &(isempty(rnap_locs) ||rnap_locs(length(rnap_locs))-RNAP_width >=0 ) ) %% make sure no rnap on the site at loading time.
			rnap_locs(length(rnap_locs)+1) = 1;
			RNAP_exitTimes(:,length(rnap_locs)) = zeros(geneLength,1);
			RNAP_RiboCoupling(length(RNAP_RiboCoupling)+1) = 0;
			loadt = t + exprnd(1/kLoading);% calculate load time for next RNAP
			Riboloadt(length(rnap_locs)) = t + exprnd(1/kRiboLoading); %calculate first ribo load time on that RNAP
			if Random_Model ==1 %do u want to use normal distribution instead of exp distribution.
				loadt= t+ normrnd(1/kRiboLoading,sqrt(1/kLoading));
				Riboloadt(length(rnap_locs)) = t + normrnd(1/kRiboLoading,sqrt(1/kRiboLoading));
			end
			Ribo_locs(length(rnap_locs),1)=0;
			deg_locs(length(rnap_locs))=0;
			riboExitTimes(:,length(rnap_locs),1)=zeros(geneLength,1);
			deg_exitTimes(:,length(rnap_locs))=zeros(geneLength,1);
			RNAlifetimes(length(rnap_locs)) = t+exprnd(mu1); %2*simtime; %calculate cootranscriptional lifetime
			RNAlifetimes_3prime(length(rnap_locs)) = 2*simtime; %t+exprnd(mu1);
			RNAlifetimes_analysis(length(rnap_locs)) = RNAlifetimes(length(rnap_locs))-t;%analysis tool
			for rs_idx =1:size(rut_sites,2)
				rho_locs(length(rnap_locs),rs_idx) =0; %initialize that this can have rut sites.
				rut_loadT(length(rnap_locs),rs_idx) = simtime +1; % initialize the loadT array but we will set the correct loadt later
			end
		end

		%%RNAP loop
		for rnap = 1:length(rnap_locs) 
			currentRNAPloc =rnap_locs(rnap);
                    
			if rnap_locs(rnap) <=geneLength
				bases_evaluated =ceil(avgSpeed1*10*dt);
				if rnap_locs(rnap)+bases_evaluated <=geneLength % add in integer(value*dt)
					if t <glutime
						tempExitTimes = t+cumsum(exprnd(specificDwelltime1(currentRNAPloc:currentRNAPloc+bases_evaluated))); %RNAP_exitTimes(currentRNAPloc,rnap)
					else
						tempExitTimes = t+cumsum(exprnd(specificDwelltime2(currentRNAPloc:currentRNAPloc+bases_evaluated)));
					end
				else
					if t< glutime
						tempExitTimes = t+cumsum(exprnd(specificDwelltime1(rnap_locs(rnap):geneLength)));
					else
						tempExitTimes = t+cumsum(exprnd(specificDwelltime2(rnap_locs(rnap):geneLength)));
					end
				end

				tempRNAP_exitTimes=  tempExitTimes((tempExitTimes>=t & tempExitTimes<=t+dt));	

				if rnap>1
					PrevRNAPloc =rnap_locs(rnap-1);
					if PrevRNAPloc == geneLength + 10
						j=1;
						while j <= size(rnap_locs(1:rnap-1)) & rnap_locs(rnap-j) == geneLength + 10 &rnap-j >1;
							j =j+1;
						end
						if j == rnap || rnap-j <1 %% if there is no rnap behind thats not terminated
							PrevRNAPloc =geneLength+1;
						else
							PrevRNAPloc = rnap_locs(rnap-j);
						end
					end
					overlap = (rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1)-PrevRNAPloc +RNAP_width; %% check for collision overlapt is +ve if collision
					
					if PrevRNAPloc >= geneLength
						overlap =0;  %if the previous rnap is done transcribing it cant collide.
                    end


					if overlap <=0 
						RNAP_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1),rnap)=tempRNAP_exitTimes;
						deg_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1),rnap) = t +exprnd(mu1*ones(length(tempRNAP_exitTimes),1));
						rnap_locs(rnap) = rnap_locs(rnap)+size(tempRNAP_exitTimes,1);
					else
						RNAP_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1-overlap),rnap)=tempRNAP_exitTimes(1:size(tempRNAP_exitTimes,1)-overlap);
						deg_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1-overlap),rnap) = t +exprnd(mu1*ones(length(tempRNAP_exitTimes)-overlap,1));
						rnap_locs(rnap) = rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-overlap;
					end
                else
					tempRNAP_exitTimes=  tempExitTimes((tempExitTimes>=t & tempExitTimes<=t+dt));
					RNAP_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1),rnap)=tempRNAP_exitTimes;
					deg_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1),rnap) = t + exprnd(mu1*ones(length(tempRNAP_exitTimes),1));
					rnap_locs(rnap) = rnap_locs(rnap)+size(tempRNAP_exitTimes,1);
				end
				if rnap_locs(rnap) >= geneLength+1 % is it done transcribing the DNA.
					if ribofollowing== 0
						condition = RNAP_exitTimes(:,rnap)>0 & deg_exitTimes(:,rnap)> t;
						deg_exitTimes(condition,rnap) = t + exprnd(mu2*ones(sum(condition),1)); 
					end
					if t < RNAlifetimes(rnap)
						RNAlifetimes(rnap) = t+ exprnd(mu2); % calculating post transcriptional degradation time.
						RNAlifetimes_analysis(rnap) = RNAlifetimes(rnap)-t; %analysis tool
						
						if rnap_locs(rnap) <= geneLength+8
							if size(Ribo_locs,1)<rnap || Ribo_locs(rnap,1) ==0 %check ribo RNAP distance.
								RiboStatistics(rnap,1) =-1;
							else
								RiboStatistics(rnap,1) = rnap_locs(rnap)-Ribo_locs(rnap,1);
							end
						else
							RiboStatistics(rnap,1) =-1;
						end
					end
					if t < RNAlifetimes_3prime(rnap)
						RNAlifetimes_3prime(rnap) = t+exprnd(mu2);
					end
				end
			end

			%%lifetime calculations for random decay model:
			if (RNAP_exitTimes(probe1,rnap) >=t &RNAP_exitTimes(probe1,rnap) <t+dt & ribofollowing ==0)
				RNAlifetimes(rnap) =t+exprnd(mu1); %calculate cootranscriptional lifetime
				RNAlifetimes_analysis(rnap) = RNAlifetimes(rnap)-t;
			elseif (RNAP_exitTimes(probe2,rnap) >=t &RNAP_exitTimes(probe2,rnap) <t+dt)
				RNAlifetimes_3prime(rnap) =t+exprnd(mu1);
			end

			%%Ribosome loading
			
			if(Riboloadt(rnap) <=t & Riboloadt(rnap) <= RNAlifetimes(rnap) & (rnap_locs(rnap) <=Ribo_width||(sum(Ribo_locs(rnap,:)>0)>0&&Ribo_locs(rnap,sum(Ribo_locs(rnap,:)>0))<= Ribo_width)) )  %% make sure no rnap on the site at loading time.
				Riboloadt(rnap) = t + exprnd(1/kRiboLoading);
			end
			
			%if(Riboloadt(rnap) <=t & Riboloadt(rnap) <= RNAlifetimes(rnap) & rnap_locs(rnap) >=Ribo_width)
			if(Riboloadt(rnap) <=t & Riboloadt(rnap) <= RNAlifetimes(rnap) & rnap_locs(rnap) >=Ribo_width & (sum(Ribo_locs(rnap,:)>0)==0||Ribo_locs(rnap,sum(Ribo_locs(rnap,:)>0))> Ribo_width)) %% check if degradation started
				Riboloadt(rnap) = t + exprnd(1/kRiboLoading);
				if loadOneRibo==1
					Riboloadt(rnap) = simtime;
				end
				if size(Ribo_locs)==0
					Ribo_locs(1,1)=0; %floor(Ribo_width/2+1);
					riboExitTimes(1,1,1) =0;
					Ribo_locs(rnap,1) = 1;
					riboExitTimes(1,rnap)= t
				elseif size(Ribo_locs,1) <rnap
					Ribo_locs(rnap,1)=1;%floor(Ribo_width/2+1);
					riboExitTimes(1,rnap,1) =t;
				else
					Ribo_locs(rnap,size(Ribo_locs(Ribo_locs(rnap,:)>=1),2)+1)=1;
					riboExitTimes(:,rnap,size(Ribo_locs(Ribo_locs(rnap,:)>=1),2)) = zeros(geneLength,1,1);
				end
            end

			%in the case of 5'->3' decay at ribo speed lets load a last ribo at the very end and use that to track degradation.
			if(t <= RNAlifetimes(rnap) & t+dt > RNAlifetimes(rnap) & ribofollowing ==2)
                if rnap_locs(rnap) <= RNAP_width ||  (sum(Ribo_locs(rnap,:)>0)>0 && Ribo_locs(rnap,sum(Ribo_locs(rnap,:)>0))<= Ribo_width)
                    RNAlifetimes(rnap) = RNAlifetimes(rnap) + dt;
                else
				    if size(deg_locs)==0
					    deg_locs(1)=0; %floor(Ribo_width/2+1);
					    deg_locs(rnap) = 1;
					    deg_exitTimes(:,1)= zeros(geneLength,1,1);
					    deg_exitTimes(1,rnap)= t;
				    else
					    deg_locs(rnap) = 1;
					    deg_exitTimes(:,rnap) = zeros(geneLength,1,1);
					    deg_exitTimes(1,rnap) =t;
                    end
                end
			end
		end

		%%rhofactor simulation.
		for RNA = 1:length(rnap_locs)

			%lets add in some cool mechanics for the PT percentage based on free RNA behind the RNAP.
			if PT_Model==2 && rnap_locs(RNA)<geneLength
				PTRNAsize = rnap_locs(RNA)- RNAP_width-rho_width;
				if sum(Ribo_locs(RNA,:)>=1) >=1
					PTRNAsize= PTRNAsize - Ribo_locs(RNA,1);
				end
				if PTRNAsize>minRholoadRNA && 100*dt*rand <= PTpercent*PTRNAsize/3075 %krutloading is calculated by using lacZ
					temp_rho_loading_loc = rnap_locs(RNA)- RNAP_width- floor(rand*PTRNAsize);
					if temp_rho_loading_loc>rho_locs(RNA,1)
						rho_locs(RNA,1)=temp_rho_loading_loc;
					end
				end
			end

			%for loop for RUT site PT cases.
			for rs_idx =1:size(rut_loadT,2) %rs_idx= rut site idx
				rut_site= rut_sites(rs_idx);
				%RUT site mechanics with Percentage Premature termination.
				if  PT_Model==1 && RNAP_exitTimes(rut_site,RNA) <t+dt && RNAP_exitTimes(rut_site,RNA) >t &100*rand <= PTpercent	
					rnap_locs(RNA)=geneLength+10;
					if(t<RNAlifetimes(RNA))
						RNAlifetimes(RNA) = t+ exprnd(muPT);
						RNAlifetimes_3prime(RNA)=t;
					end
					%rho_locs(RNA,rs_idx) =geneLength+9;
				end
				
				% RUT site mechanics with kloading
				if PT_Model==0 &&Ribo_locs(RNA,1) <= rut_site &&RNAP_exitTimes(rut_site,RNA) <t+dt && RNAP_exitTimes(rut_site,RNA) >t %%rutsite loadt calculation
					
					rut_loadT(RNA,rs_idx) = t + exprnd(1/KRutLoading);
					if Random_Model ==1 
						rut_loadT(RNA,rs_idx) = t + normrnd(1/KRutLoading,sqrt(1/KRutLoading));
					end
				end

				
				%load the rut on rutsite when loadt happens
				%%if ribo is ahead no point in having rut load so we wont include that case.
				if t>rut_loadT(RNA,rs_idx) && rho_locs(RNA,rs_idx) ==0 && Ribo_locs(RNA,1)<rut_site && rnap_locs(RNA)>rut_site && RNAP_RiboCoupling(RNA) ==0 && rnap_locs(RNA) <geneLength+1
					rho_locs(RNA,rs_idx) =rut_site;
				end

				%check if rnap already terminated
				if rnap_locs(RNA)==geneLength+10;
					rho_locs(RNA,rs_idx) =geneLength+9; % no need to have the rho continue till end lets just have it terminate
					r_loc_time(RNA,rs_idx,floor(t/dt))=0;
				end

				%%calculate rho movement.
				if(rho_locs(RNA,rs_idx) >0 && rho_locs(RNA,rs_idx) <geneLength)
                	bases_evaluated = ceil(rutSpeed*dt*10);
					if rho_locs(RNA,rs_idx)+riboSpeed*5 <=geneLength
						tempExitTimes = t+cumsum(exprnd(specificDwelltimeRho(rho_locs(RNA,rs_idx):rho_locs(RNA,rs_idx)+riboSpeed*5)));	
					else
						tempExitTimes = t+cumsum(exprnd(specificDwelltimeRho(rho_locs(RNA,rs_idx):geneLength)));	
					end
					tempRho=  tempExitTimes((tempExitTimes>=t & tempExitTimes<=t+dt));
					rho_locs(RNA,rs_idx) = rho_locs(RNA,rs_idx)+size(tempRho,1);

					r_loc_time(RNA,rs_idx,floor(t/dt))=rho_locs(RNA,rs_idx);
				end
				if rho_locs(RNA,rs_idx) >=rnap_locs(RNA)
					rnap_locs(RNA)=geneLength+10;
					condition = RNAP_exitTimes(:,RNA)>0 & deg_exitTimes(:,RNA)> t;
					deg_exitTimes(condition,RNA) = t +exprnd(muPT*ones(sum(condition),1)); 
					if(RNAlifetimes(RNA) >=t)
						RNAlifetimes(RNA) = t+ exprnd(muPT);
						RNAlifetimes_3prime(RNA)=t+exprnd(muPT);
					end
					rho_locs(RNA,rs_idx) =geneLength+9;
				end
			end
		end

		%%ribo simulation
		for RNA = 1:size(Ribo_locs,1)
			for ribo =1:size(Ribo_locs(Ribo_locs(RNA,:)>0),2)
				if(Ribo_locs(RNA,ribo) < deg_locs(RNA))
					Ribo_locs(RNA,ribo) = geneLength+10;
				end
                if Ribo_locs(RNA,ribo) <=geneLength
                	bases_evaluated = ceil(riboSpeed*10*dt);
				    if Ribo_locs(RNA,ribo)+bases_evaluated <=geneLength
					    tempExitTimes2 = t+cumsum(exprnd(RibospecificDwelltime1(Ribo_locs(RNA,ribo):Ribo_locs(RNA,ribo)+bases_evaluated)));	
				    else
					    tempExitTimes2 = t+cumsum(exprnd(RibospecificDwelltime1(Ribo_locs(RNA,ribo):geneLength)));	
				    end
				    tempRibo_exitTimes=  tempExitTimes2((tempExitTimes2>=t & tempExitTimes2<=t+dt));
				    %%add in ribosome collision
				    %%add in ribsome exit time code
				    if ribo ==1
					    if RNAP_RiboCoupling(RNA)==1 && Ribo_locs(RNA,ribo) <= geneLength-RNAP_width
						    %RNAP and RIbo are coupled so link their movements together
						    riboExitTimes(Ribo_locs(RNA,ribo):geneLength-RNAP_width,RNA,ribo) =RNAP_exitTimes(Ribo_locs(RNA,ribo)+RNAP_width:geneLength,RNA); % this is overkill u only really need to copy till current RNAPloc but this also works just computationally worse.
						    Ribo_locs(RNA,ribo) = rnap_locs(RNA)-RNAP_width;
					    elseif RNAP_RiboCoupling(RNA)==1 && Ribo_locs(RNA,ribo) > geneLength-Ribo_width && Ribo_locs(RNA,ribo) < geneLength+1
						    %compute the last bits of the ribo translating on RNAP
						    riboExitTimes(Ribo_locs(RNA,ribo):geneLength,RNA,ribo)= t+cumsum(exprnd(RibospecificDwelltime1(Ribo_locs(RNA,ribo):geneLength)));
						    Ribo_locs(RNA,ribo) = geneLength+1;
					    elseif rnap_locs(RNA) ==geneLength+10 % premature termination code for the ribosome movement
						    riboExitTimes(Ribo_locs(RNA,ribo):length(RNAP_exitTimes(:,rnap)>0),RNA,ribo)= t+cumsum(exprnd(RibospecificDwelltime1(Ribo_locs(RNA,ribo):length(RNAP_exitTimes(:,rnap)>0))));
						    idx =length((RNAP_exitTimes(RNAP_exitTimes(:,RNA)>0)))+1;
                            % could be idx = size(RNAP_exitTimes(:,RNA)>0))+1
						    riboExitTimes(idx:geneLength,RNA,ribo)= zeros(geneLength-idx+1,1);
						    Ribo_locs(RNA,ribo) = geneLength+10;
					    else
						    overlap=(Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1)-1)-rnap_locs(RNA)+RNAP_width;
						    if rnap_locs(RNA)==geneLength+1
							    overlap =0; %if rnap is done transcribing you cant overlap with the rnap
						    end
						    if overlap >0	%if collided then lets have them coupled.
							    if(rnap_locs(RNA)<=geneLength && boolRNAPRiboCoupling ==1)
								    RNAP_RiboCoupling(RNA)=1;
							    end
							    riboExitTimes(Ribo_locs(RNA,ribo):(Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1)-1-overlap),RNA,ribo)=tempRibo_exitTimes(1:size(tempRibo_exitTimes,1)-overlap);
							    Ribo_locs(RNA,ribo) = rnap_locs(RNA)-RNAP_width;
						    else
							    riboExitTimes(Ribo_locs(RNA,ribo):(Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1)-1),RNA,ribo)=tempRibo_exitTimes;
							    Ribo_locs(RNA,ribo) = Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1);
						    end
                        end
                       
				    else
					    overlap=(Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1)-1)-Ribo_locs(RNA,ribo-1)+Ribo_width;
					    if Ribo_locs(RNA,ribo-1)>=geneLength+1
						    overlap =0; %if prev ribo is done transcribing you cant overlap with the ribo
					    end
					    idx =sum(riboExitTimes(:,RNA,ribo-1)>0)+1;
					    if (Ribo_locs(RNA,ribo-1)>=geneLength+9 && Ribo_locs(RNA,ribo) >=idx)
						    Ribo_locs(RNA,ribo) = geneLength+10;
					    elseif overlap >0	%if collided then lets have them coupled.
						    riboExitTimes(Ribo_locs(RNA,ribo):(Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1)-1-overlap),RNA,ribo)=tempRibo_exitTimes(1:size(tempRibo_exitTimes,1)-overlap);
						    Ribo_locs(RNA,ribo) = Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1)-overlap;
					    else
						    riboExitTimes(Ribo_locs(RNA,ribo):(Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1)-1),RNA,ribo)=tempRibo_exitTimes;
						    Ribo_locs(RNA,ribo) = Ribo_locs(RNA,ribo)+size(tempRibo_exitTimes,1);
                        end
    
				    end
				    if (Ribo_locs(RNA,ribo) <geneLength & t>=RNAlifetimes_3prime(RNA) & ribofollowing ==0 ) %% is for random decay model.
					    Ribo_locs(RNA,ribo) =geneLength+10;
				    end
				    if (Ribo_locs(RNA,ribo) >geneLength & Ribo_locs(RNA,ribo)<geneLength+9) 
					    Protein_lifetime =[Protein_lifetime,t+exprnd(muProtein)];
					    last_idx = find(Ribo_locs(RNA,:)>0,1,'last');
					   
				    end
			    end
            end

            %lets do some degradation particle movements!
            if (t> RNAlifetimes(RNA) & deg_locs(RNA) <=geneLength & deg_locs(RNA) >0 &ribofollowing ==2)
            	bases_evaluated = ceil(degSpeed*10*dt);
			    if deg_locs(RNA)+bases_evaluated <=geneLength
				    tempExitTimes = t+cumsum(exprnd(DegspecificDwelltime1(deg_locs(RNA):deg_locs(RNA)+bases_evaluated)));	
			    else
				    tempExitTimes = t+cumsum(exprnd(DegspecificDwelltime1(deg_locs(RNA):geneLength)));	
			    end
			    tempdeg_exitTimes=  tempExitTimes((tempExitTimes>=t & tempExitTimes<=t+dt));
				
				RNAPloc = rnap_locs(RNA);
				if RNAPloc >= geneLength
					overlap=0;
				elseif RNAPloc < geneLength
					%compute the overlap to determine if there is collision
					overlap = (deg_locs(RNA)+size(tempdeg_exitTimes,1)-1)-RNAPloc + RNAP_width; %% check for collision overlap is +ve if collision
                end
				%this below segment is if you would like collision checks with Ribosome instead
				%{ 
				PrevRibo =sum(Ribo_locs(RNA,:)>0); %calculate ribo exit time
				if PrevRibo == 0
					overlap =0;  %if the there is no prev ribo then no collision can happen.
				else
					PrevRiboloc = Ribo_locs(RNA,PrevRibo); 
					overlap = (rnap_locs(RNA)+size(tempRNAP_exitTimes,1)-1)-PrevRNAPloc +RNAP_width; %% check for collision overlapt is +ve if collision
					if PrevRiboloc >= geneLength
						overlap=0;
					elseif PrevRiboloc <= geneLength
						%compute the overlap to determine if there is collision
						overlap = (deg_locs(RNA)+size(tempdeg_exitTimes,1)-1)-PrevRiboloc +Ribo_width; %% check for collision overlap is +ve if collision
	                end
	            end
	            %}
			    %%add in ribosome collision
			    %%add in ribsome exit time code
			    idx =sum(deg_exitTimes(:,RNA)>0)+1;
			    if rnap_locs(RNA) ==geneLength+10 % premature termination code for the degradation particle
			    		%if the rnap has already finished transcribing we cant really collide into anything so might as well calculate everything
					    deg_exitTimes(deg_locs(RNA):length(RNAP_exitTimes(:,RNA)>0),RNA)= t+cumsum(exprnd(DegspecificDwelltime1(deg_locs(RNA):length(RNAP_exitTimes(:,RNA)>0))));
					    idx =length((RNAP_exitTimes(RNAP_exitTimes(:,RNA)>0)))+1;
                        % could be idx = size(RNAP_exitTimes(:,RNA)>0))+1
					    deg_exitTimes(idx:geneLength,RNA)= zeros(geneLength-idx+1,1);
					    deg_locs(RNA) = geneLength+10;
			    elseif overlap >0	
				    deg_exitTimes(idx:(idx+size(tempdeg_exitTimes,1)-1-overlap),RNA)=deg_exitTimes(1:size(tempdeg_exitTimes,1)-overlap);
				    deg_locs(RNA) = deg_locs(RNA)+size(tempdeg_exitTimes,1)-overlap;
			    else
				    deg_exitTimes(deg_locs(RNA):(deg_locs(RNA)+size(tempdeg_exitTimes,1)-1),RNA)=tempdeg_exitTimes;
				    deg_locs(RNA) = deg_locs(RNA)+size(tempdeg_exitTimes,1);
                end
            end


		end
	end

	%%add in fish signals.

	%count fish
	RNAseq = zeros(geneLength,3);
	RNAseqdec = zeros(geneLength,3);
	crazyfish = zeros(geneLength,simtime);
	NETseq = zeros(geneLength,simtime);
    
  


    %remove 1st dimension from riboExitTimes to make it compatible with T_ribolocs %better implementation then squeeze()	
	squeezed_dim = size(riboExitTimes);
	squeezed_riboExitTimes = reshape(riboExitTimes(geneLength,:,:),[squeezed_dim(2:end) 1]);

	%check if ribo locs was never initialized no ribos were ever make just make a empty ribo locs this occours if rnap is never loaded
	if isempty(Ribo_locs)
		Ribo_locs =zeros(squeezed_dim(2:end));
	end

    
    for t= dt:dt:simtime
		t_index = round(t/dt);
		fish1_signal(t_index)=sum(RNAP_exitTimes(probe1,:)<t & RNAP_exitTimes(probe1,:)>0 );
		%decayed1(t_index)=0;
		fish2_signal(t_index)=sum(RNAP_exitTimes(probe2,:)<t & RNAP_exitTimes(probe2,:)>0 );
        

        %decayed2(t_index)=sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(probe2,:)'>0); %no ribo following
		%decayed2(t_index) =0;
		shared_condition =RNAlifetimes(:)<t & RNAlifetimes(:)>0;
		if ribofollowing ==1
			decayed1(t_index) =sum(shared_condition & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & deg_exitTimes(probe1,:)'<t);
			decayed2(t_index) =sum(shared_condition & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & deg_exitTimes(probe2,:)'<t );
			
        else
			decayed1(t_index) =sum(shared_condition & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & deg_exitTimes(probe1,:)'<t);
			decayed2(t_index) =sum(shared_condition & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & deg_exitTimes(probe2,:)'<t );
            %decayed1(t_index) =sum(RNAlifetimes(:)<(t-probe1/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0);
			%decayed2(t_index) =sum(RNAlifetimes(:)<(t-probe2/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 ) ;
		end
		%RNAP_exitTimes(geneLength,:)'>t
		%cootranscriptional degradation
		shared_condition=max(RNAP_exitTimes(:,:))'>t & RNAP_exitTimes(probe1,:)'<t;
		transcribing_5signal(t_index,1) = sum (shared_condition); %%total number of active RNAPs b/w probe 1 and end
		transcribing_5signal(t_index,2) = sum (shared_condition & RNAlifetimes(:)>t); %active 3' signal from transcribing rnaps
		transcribing_5signal(t_index,3) = sum (shared_condition & RNAlifetimes(:)<t); %decayed 3' signal from transcribing rnaps
		
		%premature termination
		shared_condition=RNAlifetimes(:)>t;
		prematureterm(t_index,1)= sum(shared_condition& rnap_locs(:) >= geneLength &max(RNAP_exitTimes(:,:))'<t);
		prematureterm(t_index,2)= sum(shared_condition & rnap_locs(:) ==geneLength+10 &max(RNAP_exitTimes(:,:))'<t);
		prematureterm(t_index,3)= sum(shared_condition & rnap_locs(:) < geneLength+9 &max(RNAP_exitTimes(:,:))'<t);
		
		if ribofollowing ==1
			prematureterm(t_index,1)= sum(RNAlifetimes(:)>t& rnap_locs(:) >= geneLength &max(RNAP_exitTimes(:,:))'<t);
			prematureterm(t_index,2)= sum(RNAlifetimes(:)>t & rnap_locs(:) ==geneLength+10 &max(RNAP_exitTimes(:,:))'<t);
			prematureterm(t_index,3)= sum(RNAlifetimes(:)>t & rnap_locs(:) < geneLength+9 &max(RNAP_exitTimes(:,:))'<t);
		
		end

		protein_signal(t_index) = sum(riboExitTimes(end,:,:)<t & riboExitTimes(end,:,:)>0 ,"all");
		%%protein_signal(t_index) =  protein_signal(t_index)- sum(Protein_lifetime(:)<t -sum(deg_exitTimes(:) >t));
		protein_signal(t_index) = protein_signal(t_index)- sum(Protein_lifetime(:)<t);
		%

		if (t<=glutime && (t+dt) > glutime)
			for ix =1:dx:geneLength
				shared_condition=RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)' < t;
				RNAseq(ix,1) = sum(shared_condition );
				RNAseq(ix,2) = sum(shared_condition & max(RNAP_exitTimes(:,:))'>t);
				RNAseq(ix,3) = sum(shared_condition & max(RNAP_exitTimes(:,:))'<=t);
				RNAseq(ix,4) = sum(shared_condition & max(RNAP_exitTimes(:,:))'<=t & RNAlifetimes(:)>t);
				RNAseq(ix,5) = sum(shared_condition & max(RNAP_exitTimes(:,:))'<t & RNAlifetimes(:)<t);
				if ribofollowing ==1
					shared_condition=RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & deg_exitTimes(ix,:)'<t;
					RNAseqdec(ix,1) = sum(shared_condition);
					RNAseqdec(ix,2) = sum(shared_condition & max(RNAP_exitTimes(:,:))'>t);
					RNAseqdec(ix,3) = sum(shared_condition & max(RNAP_exitTimes(:,:))'<=t);
					RNAseqdec(ix,4) = 0;
					RNAseqdec(ix,5) = sum(shared_condition & max(RNAP_exitTimes(:,:))'<=t);
				else
					shared_condition=RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & deg_exitTimes(ix,:)'<t;
					%RNAseqdec(ix,1) = sum(RNAlifetimes(:)<(t-ix/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t);
					%RNAseqdec(ix,2) = sum(RNAlifetimes(:)<(t-ix/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t& max(RNAP_exitTimes(:,:))'>t);
					%RNAseqdec(ix,3) = sum(RNAlifetimes(:)<(t-ix/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t& max(RNAP_exitTimes(:,:))'<=t);
					RNAseqdec(ix,1) = sum(shared_condition);
					RNAseqdec(ix,2) = sum(shared_condition & max(RNAP_exitTimes(:,:))'>t);
					RNAseqdec(ix,3) = sum(shared_condition & max(RNAP_exitTimes(:,:))'<=t);
					RNAseqdec(ix,4) = 0;
					RNAseqdec(ix,5) = sum(shared_condition & max(RNAP_exitTimes(:,:))'>t);
				end
			end
		end
	end
	for t= 1:simtime
		tempNETseq = sum(RNAP_exitTimes(:,:)<=t & RNAP_exitTimes(:,:)>0,1);
		tempNETseq =tempNETseq(tempNETseq>0 &tempNETseq<geneLength & max(RNAP_exitTimes(:,:),[],1)>t);
		tempNETseq = histcounts(tempNETseq,'BinMethod','integers','BinLimits',[1,geneLength]);
		NETseq(:,t)= tempNETseq;
	end
	for genex=1:geneLength
		for t =1 :simtime
			crazyfish(genex,floor(t))= sum(RNAP_exitTimes(genex,:)<t & RNAP_exitTimes(genex,:)>0 );
			if ribofollowing == 1
				crazyfish(genex,floor(t))= crazyfish(genex,floor(t)) - sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(genex,:)'<t & RNAP_exitTimes(genex,:)'>0 & deg_exitTimes(genex,:)'<t);
			else
				crazyfish(genex,floor(t))= crazyfish(genex,floor(t)) - sum(RNAlifetimes(:)<(t-genex/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(genex,:)'<t & RNAP_exitTimes(genex,:)'>0 );
			end

        	bp = genex;

		end
	end
	shared_condition = RNAP_exitTimes>0;
	SendSeq = zeros(size(RNAP_exitTimes,2),4);
	SendSeq(:,1) = sum(deg_exitTimes(:,:)<glutime & shared_condition & RNAP_exitTimes<glutime,1)';
	SendSeq(:,2) = sum(RNAP_exitTimes<glutime  & shared_condition,1)';
	SendSeq(:,3) = max(RNAP_exitTimes(:,:))'<glutime;
	SendSeq(:,4) = sum(shared_condition,1)';
	SendSeq = SendSeq(SendSeq(:,1)~=SendSeq(:,2),:);

	ribo_loadTimes=reshape(riboExitTimes(1,:,:),[squeezed_dim(2:end) 1]);
	ribo_loadTimes(:,2:end)= diff(ribo_loadTimes,1,2);
	ribo_loadTimes(:,1)= ribo_loadTimes(:,1)- RNAP_exitTimes(1,:)'; 
	statisticTerms(1) = sum(rnap_locs(:)>=geneLength);
	statisticTerms(2) = sum(Ribo_locs(:,:)==geneLength+1 ,'all');
	statisticTerms(3) = sum(rnap_locs(:)==geneLength+10);
	statisticTerms(4) = sum(Ribo_locs(:,:)>0,'all');
	fish1_signal=fish1_signal-decayed1;
	fish2_signal=fish2_signal-decayed2;


	%lets create an output struct
	output = struct;
	output.parameters = parameters;
	output.ribo_loadTimes =ribo_loadTimes(:);
	output.ribo_loadTimes2 =riboExitTimes(1,:);
	output.ribo_loadTimes2 =output.ribo_loadTimes2';
	output.fish5Signal= fish1_signal;
	output.fish3Signal = fish2_signal;
	output.proteinSignal = protein_signal;
	output.RNAlifetimes_analysis = RNAlifetimes_analysis;
	output.transcribing_5signal = transcribing_5signal;
	output.RNAP_exitTimes = RNAP_exitTimes;
	output.riboExitTimes = riboExitTimes;
	output.r_loc_time = r_loc_time;
	output.DNAStats = statisticTerms;
	output.DNAStatsLabels = ["#rna produced";"#proteins produced";"total# premature terminations";"total# of ribosomes"];
	output.FirstRiboRNAPdistance = RiboStatistics(:,1);
	output.RiboperRNAP = sum(riboExitTimes(end,:,:)>0,2);
	output.prematureterm = prematureterm;
    output.fishTime= [dt:dt:simtime]';
    output.RNAseq = RNAseq -RNAseqdec;
    output.NETseq=NETseq;
    output.Sendseq = SendSeq;
    output.RNAseqtime = glutime;
    output.crazyfish = crazyfish;
	output.RNAP_startTimes = RNAP_exitTimes(1,:);
	output.RNAP_endTimes = RNAP_exitTimes(geneLength,:);
	output.RNAlifetimes = RNAlifetimes;
	output.mu1lifetimes = RNAlifetimes(output.RNAP_endTimes(:)>RNAlifetimes(:))-output.RNAP_startTimes(output.RNAP_endTimes(:)>RNAlifetimes(:));
	output.mu2lifetimes = RNAlifetimes(output.RNAP_endTimes(:)<RNAlifetimes(:))-output.RNAP_endTimes(output.RNAP_endTimes(:)<RNAlifetimes(:));
	output.deg_exitTimes = deg_exitTimes;
    output.RNAP_dwellTimeProfile = specificDwelltime1;
    output.ribolocs = Ribo_locs;
    output.rnaplocs = rnap_locs;
    output.codelastEdited = '9th-Oct-2025';
    output.codeComment ='Changed NETseq results so that they no longer cover the RNAP footprint';

    tEnd = toc(tStart);
    output.timeElapsed = tEnd;
end
