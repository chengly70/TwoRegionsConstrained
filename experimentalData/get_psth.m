% Script to calc PSTH over all PC and OB cells

Twin=0.005; %1ms 

WinOvrLap=0;
numEvok = (2/Twin); %assuming non-overlapping windows
numSpon = (28/Twin); %discard data after 30 sec

numTot=numEvok+numSpon;

%--- output variables ----
tm_vec=(0:Twin:30-Twin)';
OB_indiv_psth=[];  %all individ OB PSTH (trial-avg)
OB_psth=[];  %pop. PSTH
n_OBtotal=0;      %total # of OB cells across all data
PC_indiv_psth=[]; %all individ PC PSTH (trial-avg)
PC_psth=[];  %pop. PSTH
n_PCtotal=0; %total # of PC cells across all data

% run processDataSelect.m
for which_set=0:1 %exp1+2 [0] or exp3+4 [1]
[stim_trials,nEpoch,nPC,nOB]=processDataSelect(Twin,WinOvrLap,which_set);

lenTrial = size(stim_trials{2}.PCcounts,1); %30 sec per trial but some extra time before start of next

%% Big array of spike counts
PCcount_allTrial_allCell = zeros(lenTrial,nEpoch-1,nPC);
for k=1:nPC
    PCcount_allTrial = zeros(lenTrial,nEpoch-1); 
    for k1=2:nEpoch
        % Read data, put into array
        [nwin,nU]=size(stim_trials{k1}.PCcounts);
        
        if (k<=nU)   % Otherwise, No spikes belonging to this unit on that trial
                     % Recorded array was set based on max UID that
                     % appeared in that trial
            if(nwin<lenTrial)
                PCcount_allTrial(1:nwin,k1-1)=stim_trials{k1}.PCcounts(1:nwin,k);
            else
                PCcount_allTrial(1:lenTrial,k1-1)=stim_trials{k1}.PCcounts(1:lenTrial,k);
            end
        end
       
    end
    PCcount_allTrial_allCell(:,:,k) = PCcount_allTrial;
end

OBcount_allTrial_allCell = zeros(lenTrial,nEpoch-1,nOB);
for k=1:nOB
    OBcount_allTrial = zeros(lenTrial,nEpoch-1); 
    for k1=2:nEpoch
        % Read data, put into array
        [nwin,nU]=size(stim_trials{k1}.OBcounts);
        
        if (k<=nU)   % Otherwise, No spikes belonging to this unit on that trial
                     % Recorded array was set based on max UID that
                     % appeared in that trial
            if(nwin<lenTrial)
                OBcount_allTrial(1:nwin,k1-1)=stim_trials{k1}.OBcounts(1:nwin,k);
            else
                OBcount_allTrial(1:lenTrial,k1-1)=stim_trials{k1}.OBcounts(1:lenTrial,k);
            end
        end
       
    end
    OBcount_allTrial_allCell(:,:,k) = OBcount_allTrial;
end

%% PC
PCmeankeep = [];
PCmnStim1 = []; %first stimulus (and 28 sec of spont after)
PCmnStim2 = []; %second stimulus (and 28 sec of spont after)
for k=1:nPC
    PCcount_allTrial = PCcount_allTrial_allCell(:,:,k);
    
    PCmeankeep = [PCmeankeep; mean(PCcount_allTrial')];
    % Also look at first half of trials, 2nd half
    PCmnStim1 = [PCmnStim1; mean(PCcount_allTrial(:,1:10)')];
    PCmnStim2 = [PCmnStim2; mean(PCcount_allTrial(:,11:20)')];
end

%% OB
OBmeankeep = [];
OBmnStim1 = []; %first stimulus (and 28 sec of spont after)
OBmnStim2 = []; %second stimulus (and 28 sec of spont after)
for k=1:nOB
    OBcount_allTrial = OBcount_allTrial_allCell(:,:,k);
    
    OBmeankeep = [OBmeankeep; mean(OBcount_allTrial')];
    % Also look at first half of trials, 2nd half
    OBmnStim1 = [OBmnStim1; mean(OBcount_allTrial(:,1:10)')];
    OBmnStim2 = [OBmnStim2; mean(OBcount_allTrial(:,11:20)')];
end

%% Exclude firing rates that are identically 0
i_PCs=(sum(PCmeankeep')~=0);     %also used in loop below for Cov/corr
PCmeankeep=PCmeankeep(i_PCs',:);
PCmnStim1=PCmnStim1(i_PCs',:);
PCmnStim2=PCmnStim2(i_PCs',:);
nPC=sum(i_PCs); %update size of PC
i_OBs=(sum(OBmeankeep')~=0);    %also used in loop below for Cov/corr
OBmeankeep=OBmeankeep(i_OBs',:);
OBmnStim1=OBmnStim1(i_OBs',:);
OBmnStim2=OBmnStim2(i_OBs',:);
nOB=sum(i_OBs); %update size of OB

%update total number of OB/PC cells
n_PCtotal=n_PCtotal+nPC;
n_OBtotal=n_OBtotal+nOB;
% Calculate firing rate stats (in time): psth, etc
PC_indiv_psth=[PC_indiv_psth; PCmeankeep(:,1:numTot)];

OB_indiv_psth=[OB_indiv_psth; OBmeankeep(:,1:numTot)];


end

PC_psth=sum(PC_indiv_psth)'./(n_PCtotal*Twin);
OB_psth=sum(OB_indiv_psth)'./(n_OBtotal*Twin);

figure
plot(tm_vec,PC_psth,'color',[0 0.5 0],'LineWidth',4)
%shadedErrorBar(tm_vec,PC_psth,std(PC_indiv_psth),{'color',[0 .5 0]},1)
box off
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('PC Firing Rate (Hz)')

figure
plot(tm_vec,OB_psth,'b','LineWidth',4)
%shadedErrorBar(tm_vec,OB_psth,std(OB_indiv_psth),'b',1)
box off
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('OB Firing Rate (Hz)')

