function output= albur_TASEP_function(input_parameters)
	%%Parameter array
	%parameter=[avgSpeed1,avgSpeed2,kLoading,simTime,gluTime,riboSpeed,KRiboLoading,KRutLoading,mu1,mu2,muPT,muprotein,ribofollowing,RNAP dwell model ,	random statistics model ,	boolRNAPRiboCoupling , bool PT%]
	%				1 		  2 		3  		4 		5 		6  		7  				8  	     9   10  11       12        13			,		14 ,				15 			,		16 				 ,   17 ]
	%parameters =[21;9.3;1/20;2000;75;21;0.2;0;600;120;120;50;1;0;0;0;1]; 
	
    tStart = tic;

	parameters =struct;parameters.avgSpeed1 =19;parameters.avgSpeed2 = 11;parameters.ribospeed= 19;
	parameters.kLoading =1/20;parameters.kRiboLoading = 0.1;parameters.KRutLoading =0;
	parameters.simtime=2000;parameters.glutime=1600; 
	parameters.mu1=600;parameters.mu2=120;parameters.muPT=120;parameters.muProtein=6000;
	parameters.ribofollowing = 2;parameters.RNAP_Profile =0;parameters.Random_Model =0;
	parameters.boolRNAPRiboCoupling = 1;parameters.PT_Model = 1; 
	parameters.GstarRatio = 0.2; parameters.GstarAdditionTime = 1500; parameters.GstarRemovalTime = 1900; 
	parameters.geneLength=3075;%DNA length
	if nargin >=1 %number of arguments in
		parameterlabels= fieldnames(input_parameters);
		for i = 1:numel(parameterlabels)
				parameters.(parameterlabels{i})= input_parameters.(parameterlabels{i});
		end
	end
	parameters.rutSpeed=5*parameters.ribospeed;


	avgSpeed1=parameters.avgSpeed1;avgSpeed2=parameters.avgSpeed2;riboSpeed = parameters.ribospeed;
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
    GstarRatio = parameters.GstarRatio; GstarAdditionTime = parameters.GstarAdditionTime; GstarRemovalTime = parameters.GstarRemovalTime;

	supercoilingmodel = 0;
	muSuperCoil=500;
	supercoilinglifetime = [glutime + exprnd(muSuperCoil)];


	% statisticTerms[1] =#rna produced statisticTerms[2] =#proteins produced
	statisticTerms =-1*ones(4,1);

	%RiboStatistics[1] = ribo-rnap distances RiboStatistics[2] = ribo per rnap 
	RiboStatistics = -1*ones(2,2);


	% Input parameters for ribosome initiation and elongation
	riboExitTimes =zeros(geneLength/dx,2,2);; % lets make it a 3d array such as [bp,RNAP#,Ribo#]
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

	tempExitTimes =[];
	%RNAlifetimes=exprnd(mu1);
	RNAlifetimes=[];
	RNAlifetimes_3prime=[];	 % remove all mentions of this to go back to 3' and 5' decay together
	RNAlifetimes_analysis =[];
	avgDwelltime1 = dx/avgSpeed1; %sec per nucleotide
	avgDwelltime2 = dx/avgSpeed2; %sec per nucleotide
	riboavgDwelltime = dx/riboSpeed;
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
	RNAP_exitTimes = zeros(geneLength/dx,1);
	tempExitTimes = cumsum(exprnd(specificDwelltime1));

    % G* things
    RNA_Gstar = zeros(1,4); 
    Gstar_Signal = zeros(simtime/dt,1);
    decayed_Gstar = zeros(simtime/dt,1);
    Gstar_fish_signals =zeros(simtime/dt,2);
	decayed_Gstar_fish_signals =zeros(simtime/dt,2);

	%%time to prep that protein
	Protein_lifetime= [];
	protein_signal=zeros(simtime/dt,1);
	isRiboReal = [];
	ribo_loadTimes = -ones(2,2);



	%Simulation Loop
	for t = 0:dt:simtime
		%%RNAP loading
		%%lets test some random bs
		
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
			riboExitTimes(:,length(rnap_locs),1)=zeros(geneLength,1);
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
				if rnap_locs(rnap)+avgSpeed1*3 <=geneLength % add in integer(value*dt)
					if t <glutime
						tempExitTimes = t+cumsum(exprnd(specificDwelltime1(rnap_locs(rnap):rnap_locs(rnap)+avgSpeed1*3))); %RNAP_exitTimes(currentRNAPloc,rnap)
					else
						tempExitTimes = t+cumsum(exprnd(specificDwelltime2(rnap_locs(rnap):rnap_locs(rnap)+avgSpeed1*3)));
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

                    if t>GstarAdditionTime && t<GstarRemovalTime
                        bpmoved = length(tempRNAP_exitTimes);
                        if overlap > 0
                            bpmoved = bpmoved - overlap;
                        end
                        for i=1:1:bpmoved
                            isBaseGstar = rand < 0.25*GstarRatio; % 0.25 corresponds to the prob of the base being G
                            if isBaseGstar
                                if size(RNA_Gstar,1) < rnap || RNA_Gstar(rnap,1) == 0 % if row corresponding to current rnap doesn't exist, or if this is the first G*
                                    RNA_Gstar(rnap,2) = tempRNAP_exitTimes(i); % record time if this is first G*
                                    RNA_Gstar(rnap,3) = currentRNAPloc + i - 1; % record first G* position
                                end
                                RNA_Gstar(rnap,1) = RNA_Gstar(rnap,1) + 1; % add one G* to current count
                                RNA_Gstar(rnap,4) = currentRNAPloc + i - 1; % record current G* position as last position % currentRNAPloc = bp currently being transcribed, use PrevRNAPloc to get bp that has finished transcribing
                            end
                        end
                    end

					if overlap <=0 
						RNAP_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1),rnap)=tempRNAP_exitTimes;
						rnap_locs(rnap) = rnap_locs(rnap)+size(tempRNAP_exitTimes,1);
					else
						RNAP_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1-overlap),rnap)=tempRNAP_exitTimes(1:size(tempRNAP_exitTimes,1)-overlap);
						rnap_locs(rnap) = rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-overlap;
					end
                else
					tempRNAP_exitTimes=  tempExitTimes((tempExitTimes>=t & tempExitTimes<=t+dt));
					RNAP_exitTimes(rnap_locs(rnap):(rnap_locs(rnap)+size(tempRNAP_exitTimes,1)-1),rnap)=tempRNAP_exitTimes;
					rnap_locs(rnap) = rnap_locs(rnap)+size(tempRNAP_exitTimes,1);
				end
				if rnap_locs(rnap) >= geneLength+1 % is it done transcribing the DNA.
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
			%%lets test some random bs
			
			if(Riboloadt(rnap) <=t & Riboloadt(rnap) <= RNAlifetimes(rnap) & (rnap_locs(rnap) <=Ribo_width||(sum(Ribo_locs(rnap,:)>0)>0&&Ribo_locs(rnap,sum(Ribo_locs(rnap,:)>0))<= Ribo_width)) )  %% make sure no rnap on the site at loading time.
				Riboloadt(rnap) = t + exprnd(1/kRiboLoading);
			end
			
			%if(Riboloadt(rnap) <=t & Riboloadt(rnap) <= RNAlifetimes(rnap) & rnap_locs(rnap) >=Ribo_width)
			if(Riboloadt(rnap) <=t & Riboloadt(rnap) <= RNAlifetimes(rnap) & rnap_locs(rnap) >=Ribo_width & (sum(Ribo_locs(rnap,:)>0)==0||Ribo_locs(rnap,sum(Ribo_locs(rnap,:)>0))> Ribo_width)) %% check if degradation started
				Riboloadt(rnap) = t + exprnd(1/kRiboLoading);
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
				    if size(Ribo_locs)==0
					    Ribo_locs(1,1)=0; %floor(Ribo_width/2+1);
					    riboExitTimes(1,1,1) =0;
					    Ribo_locs(rnap,1) = 1;
					    riboExitTimes(1,rnap)= t;
				    elseif size(Ribo_locs,1) <rnap
					    Ribo_locs(rnap,1)=1;%floor(Ribo_width/2+1);
					    riboExitTimes(1,rnap,1) =t+dt;
				    else
					    Ribo_locs(rnap,size(Ribo_locs(Ribo_locs(rnap,:)>=1),2)+1)=1;
					    riboExitTimes(:,rnap,size(Ribo_locs(Ribo_locs(rnap,:)>=1),2)) = zeros(geneLength,1,1);
					    riboExitTimes(1,rnap,size(Ribo_locs(Ribo_locs(rnap,:)>=1),2)) =t+dt;
                    end
                end
			end
		end

		%%rhofactor simulation.
		for RNA = 1:length(rnap_locs)
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

				%lets add in some cool mechanics for the PT percentage based on free RNA behind the RNAP.
				if PT_Model==2 && rnap_locs(RNA)<geneLength
					PTRNAsize = rnap_locs(RNA)- RNAP_width-rho_width;
					if sum(Ribo_locs(RNA,:)>=1) >=1
						PTRNAsize= PTRNAsize - Ribo_locs(RNA,1);
					end
					if PTRNAsize>minRholoadRNA && 100*dt*rand <= PTpercent*PTRNAsize/geneLength
						temp_rho_loading_loc = rnap_locs(RNA)- RNAP_width- floor(rand*PTRNAsize);
						if temp_rho_loading_loc>rho_locs(RNA,1)
							rho_locs(RNA,1)=temp_rho_loading_loc;
						end
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
                if Ribo_locs(RNA,ribo) <=geneLength
				    if Ribo_locs(RNA,ribo)+riboSpeed*3 <=geneLength
					    tempExitTimes2 = t+cumsum(exprnd(RibospecificDwelltime1(Ribo_locs(RNA,ribo):Ribo_locs(RNA,ribo)+riboSpeed*3)));	
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
					    if t >= RNAlifetimes(RNA) && ribo == last_idx % think about this implementation later
						    isRiboReal(end+1) = -1; % maybe use size(Ribo_locs(Ribo_locs(rnap,:)>=1),2)+1 as the check.
					    else
						    isRiboReal(end+1) = 1;
					    end
				    end
			    end
            end
		end
	end

	%%add in fish signals.

	%create lastriboexitTime
	lastRiboExitTime=zeros(geneLength,size(RNAP_exitTimes,2));
	for RNA =1:size(Ribo_locs,1)
		last_idx = find(Ribo_locs(RNA,:)>0,1,'last');
        % lastRiboExitTime(:,RNA)= riboExitTimes(:,RNA,last_idx);
		if length(last_idx ==1)   % necessary condition - if no ribos load at all, 'find' will not return anything --> last_idx is an empty variable
			lastRiboExitTime(:,RNA)= riboExitTimes(:,RNA,last_idx);
		end
	end

	%count fish
	RNAseq = zeros(geneLength,3);
	RNAseqdec = zeros(geneLength,3);
	crazyfish = zeros(geneLength,simtime);
	NETseq = zeros(geneLength,simtime);
    
    % make sure length of rnap_locs matches the size of RNA_Gstar
    if length(rnap_locs) > size(RNA_Gstar,1)
        RNA_Gstar(length(rnap_locs),1) = 0;
    end

    Gstar_DecayTime = zeros(size(RNA_Gstar,1),1);
    for i=1:1:size(RNA_Gstar,1)
        if (RNA_Gstar(i,4)>0)
            Gstar_DecayTime(i) = lastRiboExitTime(RNA_Gstar(i,4),i); % RNA_Gstar(i,4) not all positive ints, ones in beg. are 0
            if (Gstar_DecayTime(i)==0)
                % Gstar_DecayTime(i) = RNA_Gstar(i,2);
                Gstar_DecayTime(i) = simtime+1;
            end
        end
    end


    %remove 1st dimension from riboExitTimes to make it compatible with T_ribolocs %better implementation then squeeze()	
	squeezed_dim = size(riboExitTimes);
	squeezed_riboExitTimes = reshape(riboExitTimes(geneLength,:,:),[squeezed_dim(2:end) 1]);

	%check if ribo locs was never initialized no ribos were ever make just make a empty ribo locs this occours if rnap is never loaded
	if isempty(Ribo_locs)
		Ribo_locs =zeros(squeezed_dim(2:end)); T_ribolocs = zeros(squeezed_dim(2:end));T_ribolocs(:,1)= 1;
	end

    %%lets count how much protein we got when
	T_ribolocs =(Ribo_locs >0);

	if ribofollowing ==2
		last_idx = sum(T_ribolocs,2);
		if sum(T_ribolocs,2) ==0 %in case no rnaps load
			last_idx =1;
		end
		for ii=1: length(last_idx)
			if last_idx(ii)>=1
				T_ribolocs(ii,last_idx(ii)) = 0;
			end
		end
	end
	gstar_5_breakdown= zeros(simtime/dt,3);
    gstar_3_breakdown= zeros(simtime/dt,3);
    gstar_bp1_breakdown= zeros(simtime/dt,3);
    crazyGstar = zeros(simtime, geneLength);
    for t= dt:dt:simtime
		t_index = round(t/dt);
		fish1_signal(t_index)=sum(RNAP_exitTimes(probe1,:)<t & RNAP_exitTimes(probe1,:)>0 );
		%decayed1(t_index)=0;
		fish2_signal(t_index)=sum(RNAP_exitTimes(probe2,:)<t & RNAP_exitTimes(probe2,:)>0 );
        
        % the number of RNA that have received their first G*
        Gstar_Signal(t_index) = sum(RNA_Gstar(:,2)<t & RNA_Gstar(:,2)>0); 
        decayed_Gstar(t_index) = sum(Gstar_DecayTime(:)<t & Gstar_DecayTime(:)>0);
		Gstar_fish_signals(t_index,1) = sum(Gstar_DecayTime(:)'>t & RNA_Gstar(:,1)'>0 & RNA_Gstar(:,2)'<t & RNAP_exitTimes(probe1,:)<t & RNAP_exitTimes(probe1,:)>0);
        Gstar_fish_signals(t_index,2) = sum(Gstar_DecayTime(:)'>t & RNA_Gstar(:,1)'>0 & RNA_Gstar(:,2)'<t & RNAP_exitTimes(probe2,:)<t & RNAP_exitTimes(probe2,:)>0);
        
        %Gstar breakdowns:
        % % no RNAlifetimes check
        % gstar_bp1_breakdown(t_index,1) =sum(RNAP_exitTimes(1,:)'<t & RNAP_exitTimes(1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(1,:)'>=t);
        % gstar_5_breakdown(t_index,1) =sum(RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(probe1,:)'>=t);
	  	% gstar_3_breakdown(t_index,1) =sum(RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(probe2,:)'>=t);
        % gstar_bp1_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(1,:)'>=t);
        % gstar_5_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(probe1,:)'>=t);
	    % gstar_3_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(probe2,:)'>=t);
        % gstar_bp1_breakdown(t_index,3) =sum(RNAP_exitTimes(1,:)'<t & RNAP_exitTimes(1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(1,:)'>=t);
	    % gstar_5_breakdown(t_index,3) =sum(RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(probe1,:)'>=t);
		% gstar_3_breakdown(t_index,3) =sum(RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(probe2,:)'>=t);
       
        % % RNAlifetimes check
        % gstar_bp1_breakdown(t_index,1) =sum(RNAP_exitTimes(1,:)'<t & RNAP_exitTimes(1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(1,:)'>=t | RNAlifetimes(:)>=t) );
        % gstar_5_breakdown(t_index,1) =sum(RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe1,:)'>=t | RNAlifetimes(:)>=t) );
  	    % gstar_3_breakdown(t_index,1) =sum(RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe2,:)'>=t | RNAlifetimes(:)>=t) );
        % gstar_bp1_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(1,:)'>=t | RNAlifetimes(:)>=t) );
        % gstar_5_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe1,:)'>=t | RNAlifetimes(:)>=t) );
        % gstar_3_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe2,:)'>=t | RNAlifetimes(:)>=t) );
        % gstar_bp1_breakdown(t_index,3) =sum(RNAP_exitTimes(1,:)'<t & RNAP_exitTimes(1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(1,:)'>=t | RNAlifetimes(:)>=t) );
        % gstar_5_breakdown(t_index,3) =sum(RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe1,:)'>=t | RNAlifetimes(:)>=t) );
	    % gstar_3_breakdown(t_index,3) =sum(RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe2,:)'>=t | RNAlifetimes(:)>=t) );

        %	( input.lastRiboExitTime(probe1,:)'==0  |  input.lastRiboExitTime(probe1,:)'>=t  | input.RNAlifetimes(:)>=t )
        % another fix 7/29/25
        gstar_bp1_breakdown(t_index,1) =sum(RNAP_exitTimes(1,:)'<t & RNAP_exitTimes(1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(1,:)'==0  |  lastRiboExitTime(1,:)'>=t | RNAlifetimes(:)>=t) );
        gstar_5_breakdown(t_index,1) =sum(RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe1,:)'==0  |  lastRiboExitTime(probe1,:)'>=t | RNAlifetimes(:)>=t) );
  	    gstar_3_breakdown(t_index,1) =sum(RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'>t & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe2,:)'==0  |  lastRiboExitTime(probe2,:)'>=t | RNAlifetimes(:)>=t) );
        gstar_bp1_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(1,:)'==0  |  lastRiboExitTime(1,:)'>=t | RNAlifetimes(:)>=t) );
        gstar_5_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe1,:)'==0  |  lastRiboExitTime(probe1,:)'>=t | RNAlifetimes(:)>=t) );
        gstar_3_breakdown(t_index,2) =sum(RNAP_exitTimes(end,:)'<=t & RNAP_exitTimes(end,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe2,:)'==0  |  lastRiboExitTime(probe2,:)'>=t | RNAlifetimes(:)>=t) );
        gstar_bp1_breakdown(t_index,3) =sum(RNAP_exitTimes(1,:)'<t & RNAP_exitTimes(1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(1,:)'==0  |  lastRiboExitTime(1,:)'>=t | RNAlifetimes(:)>=t) );
        gstar_5_breakdown(t_index,3) =sum(RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe1,:)'==0  |  lastRiboExitTime(probe1,:)'>=t | RNAlifetimes(:)>=t) );
	    gstar_3_breakdown(t_index,3) =sum(RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & max(RNAP_exitTimes(:,:),[],1)'<=t & rnap_locs(:)>geneLength+1 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & ( lastRiboExitTime(probe2,:)'==0  |  lastRiboExitTime(probe2,:)'>=t | RNAlifetimes(:)>=t) );


        %decayed2(t_index)=sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(probe2,:)'>0); %no ribo following
		%decayed2(t_index) =0;
		if ribofollowing ==1
			decayed1(t_index) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & lastRiboExitTime(probe1,:)'<t);
			decayed2(t_index) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & lastRiboExitTime(probe2,:)'<t );
			% decayed_Gstar_fish_signals(t_index,1) =sum(Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & lastRiboExitTime(probe1,:)'<t);
            % decayed_Gstar_fish_signals(t_index,2) =sum(Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & lastRiboExitTime(probe2,:)'<t);
		    decayed_Gstar_fish_signals(t_index,1) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & lastRiboExitTime(probe1,:)'<t & lastRiboExitTime(probe1,:)'>0);
            decayed_Gstar_fish_signals(t_index,2) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & lastRiboExitTime(probe2,:)'<t & lastRiboExitTime(probe2,:)'>0);
        elseif ribofollowing ==2
			decayed1(t_index) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & lastRiboExitTime(probe1,:)'<t);
			decayed2(t_index) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & lastRiboExitTime(probe2,:)'<t ) ;
			% decayed_Gstar_fish_signals(t_index,1) =sum(Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & lastRiboExitTime(probe1,:)'<t);
            % decayed_Gstar_fish_signals(t_index,2) =sum(Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & lastRiboExitTime(probe2,:)'<t);
			decayed_Gstar_fish_signals(t_index,1) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0 & lastRiboExitTime(probe1,:)'<t & lastRiboExitTime(probe1,:)'>0);
            decayed_Gstar_fish_signals(t_index,2) =sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,1)>0 & RNA_Gstar(:,2)<t & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 & lastRiboExitTime(probe2,:)'<t & lastRiboExitTime(probe2,:)'>0);
            
            %decayed1(t_index) =sum(RNAlifetimes(:)<(t-probe1/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(probe1,:)'<t & RNAP_exitTimes(probe1,:)'>0);
			%decayed2(t_index) =sum(RNAlifetimes(:)<(t-probe2/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(probe2,:)'<t & RNAP_exitTimes(probe2,:)'>0 ) ;
		else
			decayed2(t_index)=sum(RNAlifetimes_3prime(:)<t & RNAlifetimes_3prime(:)>0 & RNAP_exitTimes(probe2,:)'>0 & RNAP_exitTimes(probe2,:)'<t); %no ribo following differe lifetime for 3'
			decayed1(t_index)=sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(probe1,:)'>0 & RNAP_exitTimes(probe1,:)'<t); %no ribo following
		end
		%RNAP_exitTimes(geneLength,:)'>t
		%cootranscriptional degradation
		transcribing_5signal(t_index,1) = sum (max(RNAP_exitTimes(:,:))'>t & RNAP_exitTimes(probe1,:)'<t); %%total number of active RNAPs b/w probe 1 and end
		transcribing_5signal(t_index,2) = sum (max(RNAP_exitTimes(:,:))'>t & RNAP_exitTimes(probe1,:)'<t & RNAlifetimes(:)>t); %active 3' signal from transcribing rnaps
		transcribing_5signal(t_index,3) = sum (max(RNAP_exitTimes(:,:))'>t & RNAP_exitTimes(probe1,:)'<t & RNAlifetimes(:)<t); %decayed 3' signal from transcribing rnaps
		
		%premature termination
		prematureterm(t_index,1)= sum(RNAlifetimes(:)>t& rnap_locs(:) >= geneLength &max(RNAP_exitTimes(:,:))'<t);
		prematureterm(t_index,2)= sum(RNAlifetimes(:)>t & rnap_locs(:) ==geneLength+10 &max(RNAP_exitTimes(:,:))'<t);
		prematureterm(t_index,3)= sum(RNAlifetimes(:)>t & rnap_locs(:) < geneLength+9 &max(RNAP_exitTimes(:,:))'<t);
		
		if ribofollowing ==1
			prematureterm(t_index,1)= sum(RNAlifetimes(:)>t& rnap_locs(:) >= geneLength &max(RNAP_exitTimes(:,:))'<t);
			prematureterm(t_index,2)= sum(RNAlifetimes(:)>t & rnap_locs(:) ==geneLength+10 &max(RNAP_exitTimes(:,:))'<t);
			prematureterm(t_index,3)= sum(RNAlifetimes(:)>t & rnap_locs(:) < geneLength+9 &max(RNAP_exitTimes(:,:))'<t);
		
		end

		protein_signal(t_index) = sum(squeezed_riboExitTimes(:,:)<t & squeezed_riboExitTimes(:,:)>0 & T_ribolocs(:,:)==1,"all");
		%%protein_signal(t_index) =  protein_signal(t_index)- sum(Protein_lifetime(:)<t -sum(lastRiboExitTime(:) >t));
		protein_signal(t_index) = protein_signal(t_index)- sum(Protein_lifetime(:)<t & isRiboReal(:) ==1,'all');
		%

		if (t<=glutime && (t+dt) > glutime)
			for ix =1:dx:geneLength
				RNAseq(ix,1) = sum(RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)' <t);
				RNAseq(ix,2) = sum(RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)' <t & max(RNAP_exitTimes(:,:))'>t);
				RNAseq(ix,3) = sum(RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)' <t & max(RNAP_exitTimes(:,:))'<=t);
				RNAseq(ix,4) = sum(RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)' <t & max(RNAP_exitTimes(:,:))'<=t & RNAlifetimes(:)>t);
				RNAseq(ix,5) = sum(RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)' <t & max(RNAP_exitTimes(:,:))'<t & RNAlifetimes(:)<t);
				if ribofollowing ==1
					RNAseqdec(ix,1) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t);
					RNAseqdec(ix,2) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t & max(RNAP_exitTimes(:,:))'>t);
					RNAseqdec(ix,3) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t & max(RNAP_exitTimes(:,:))'<=t);
					RNAseqdec(ix,4) = 0;
					RNAseqdec(ix,5) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t & max(RNAP_exitTimes(:,:))'<=t);
				elseif ribofollowing ==2
					%RNAseqdec(ix,1) = sum(RNAlifetimes(:)<(t-ix/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t);
					%RNAseqdec(ix,2) = sum(RNAlifetimes(:)<(t-ix/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t& max(RNAP_exitTimes(:,:))'>t);
					%RNAseqdec(ix,3) = sum(RNAlifetimes(:)<(t-ix/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t& max(RNAP_exitTimes(:,:))'<=t);
					RNAseqdec(ix,1) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t);
					RNAseqdec(ix,2) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t & max(RNAP_exitTimes(:,:))'>t);
					RNAseqdec(ix,3) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t & max(RNAP_exitTimes(:,:))'<=t);
					RNAseqdec(ix,4) = 0;
					RNAseqdec(ix,5) = sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(ix,:)'>0 & RNAP_exitTimes(ix,:)'<t & lastRiboExitTime(ix,:)'<t & max(RNAP_exitTimes(:,:))'>t);
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
				crazyfish(genex,floor(t))= crazyfish(genex,floor(t)) - sum(RNAlifetimes(:)<t & RNAlifetimes(:)>0 & RNAP_exitTimes(genex,:)'<t & RNAP_exitTimes(genex,:)'>0 & lastRiboExitTime(genex,:)'<t);
			elseif ribofollowing ==2
				crazyfish(genex,floor(t))= crazyfish(genex,floor(t)) - sum(RNAlifetimes(:)<(t-genex/riboSpeed) & RNAlifetimes(:)>0 & RNAP_exitTimes(genex,:)'<t & RNAP_exitTimes(genex,:)'>0 );
			end
			%makes RNAP take the size of RNAP width in NETSEQ but that should not be the case
			%{
			if (genex <geneLength- RNAP_width+1)
				NETseq(genex,t) = NETseq(genex,t) +sum(NETseq(genex+1:genex+RNAP_width-1,t));
			elseif(genex <geneLength)
				NETseq(genex,t) = NETseq(genex,t) +sum(NETseq(genex+1:end,t));
            end
            %} 

        	bp = genex;
            crazyGstar(t, bp) = sum( RNAP_exitTimes(bp,:)' <t & RNAP_exitTimes(bp,:)'>0 & Gstar_DecayTime(:)>t & RNA_Gstar(:,2)<t & lastRiboExitTime(bp,:)'>=t ) ;

		end
	end
	SendSeq = zeros(size(RNAP_exitTimes,2),3);
	SendSeq(:,1) = sum(lastRiboExitTime(:,:)<glutime & RNAP_exitTimes>0 & RNAP_exitTimes<glutime,1)';
	SendSeq(:,2) = sum(RNAP_exitTimes<glutime  & RNAP_exitTimes>0,1)';
	SendSeq(:,3) = max(RNAP_exitTimes(:,:))'<glutime;
	SendSeq(:,4) = sum(RNAP_exitTimes>0,1)';

	ribo_loadTimes=reshape(riboExitTimes(1,:,:),[squeezed_dim(2:end) 1]);
	ribo_loadTimes(:,2:end)= diff(ribo_loadTimes,1,2);
	ribo_loadTimes(:,1)= ribo_loadTimes(:,1)- RNAP_exitTimes(1,:)';
	statisticTerms(1) = sum(rnap_locs(:)>=geneLength);
	statisticTerms(2) = sum(Ribo_locs(:,:)==geneLength+1 ,'all');
	statisticTerms(3) = sum(rnap_locs(:)==geneLength+10);
	statisticTerms(4) = sum(Ribo_locs(:,:)>0,'all');
	fish1_signal=fish1_signal-decayed1;
	fish2_signal=fish2_signal-decayed2;
    Gstar_Signal = Gstar_Signal - decayed_Gstar;
    Gstar_fish_signals = Gstar_fish_signals - decayed_Gstar_fish_signals;


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
	output.RiboperRNAP = sum(T_ribolocs(:,:)>0,2);
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
	output.RNA_Gstar = RNA_Gstar;
    output.Gstar_Signal = Gstar_Signal;
    output.Gstar_fishSignals = Gstar_fish_signals;
    output.gstar_5_breakdown= gstar_5_breakdown;
    output.gstar_3_breakdown= gstar_3_breakdown;
    output.gstar_bp1_breakdown= gstar_bp1_breakdown;
    output.lastRiboExistTime = lastRiboExitTime;
    output.crazyGstar = crazyGstar;
    output.RNAP_dwellTimeProfile = specificDwelltime1;
    output.ribolocs = Ribo_locs;
    output.rnaplocs = rnap_locs;
    output.codelastEdited = '9th-Oct-2025';
    output.codeComment ='Changed NETseq results so that they no longer cover the RNAP footprint';

    tEnd = toc(tStart);
    output.timeElapsed = tEnd;
end

