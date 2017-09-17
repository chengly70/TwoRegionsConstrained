% Script to calc firing rates; uses processDataSelect.m

which_set=input('Enter 0 (exp1+2) or 1 (exp3+4): '); 

Twin=2; %(sec) DO NOT CHANGE!
WinOvrLap=0; %0=disjoint windows, 1=overlapping windows
if(Twin==2 || WinOvrLap==0) %if Twin=2, no choice, largest length of evoked window
    WinOvrLap=0;
    numEvok = (2/Twin); %assuming non-overlapping windows
    numSpon = (28/Twin); %discard data after 30 sec
else
    numEvok = 2*(2/Twin)-1; 
    numSpon = 2*(28/Twin); %first window can technically be in 2sec 'evoked' time
end
numTot=numEvok+numSpon;

% run processDataSelect.m
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
%PCstd = []; %std across trial
for k=1:nPC
    
    PCcount_allTrial = PCcount_allTrial_allCell(:,:,k);
    PCmeankeep = [PCmeankeep; mean(PCcount_allTrial')];
    %PCstd  = [PCstd; std(PCcount_allTrial')];
    
    % Also look at first half of trials, 2nd half
    PCmnStim1 = [PCmnStim1; mean(PCcount_allTrial(:,1:10)')];
    PCmnStim2 = [PCmnStim2; mean(PCcount_allTrial(:,11:20)')];
end

%% OB
OBmeankeep = [];
OBmnStim1 = []; %first stimulus (and 28 sec of spont after)
OBmnStim2 = []; %second stimulus (and 28 sec of spont after)
%OBstd = []; %std across trial

for k=1:nOB
    
    OBcount_allTrial = OBcount_allTrial_allCell(:,:,k);
    OBmeankeep = [OBmeankeep; mean(OBcount_allTrial')];
    %OBstd  = [OBstd; std(OBcount_allTrial')];
    
    % Also look at first half of trials, 2nd half
    OBmnStim1 = [OBmnStim1; mean(OBcount_allTrial(:,1:10)')];
    OBmnStim2 = [OBmnStim2; mean(OBcount_allTrial(:,11:20)')];
end

%% Exclude firing rates that are identically 0
i_PCs=(sum(PCmeankeep')~=0);
PCmeankeep=PCmeankeep(i_PCs',:);
PCmnStim1=PCmnStim1(i_PCs',:);
PCmnStim2=PCmnStim2(i_PCs',:);
nPC=sum(i_PCs); %update size of PC
i_OBs=(sum(OBmeankeep')~=0);
OBmeankeep=OBmeankeep(i_OBs',:);
OBmnStim1=OBmnStim1(i_OBs',:);
OBmnStim2=OBmnStim2(i_OBs',:);
nOB=sum(i_OBs); %update size of OB

% Translate into firing rates
PC_sponBoth=sum(PCmeankeep(:,numEvok+1:numTot),2)./(numSpon*Twin);
PC_odBoth=sum(PCmeankeep(:,1:numEvok),2)./(numEvok*Twin);
OB_sponBoth=sum(OBmeankeep(:,numEvok+1:numTot),2)./(numSpon*Twin);
OB_odBoth=sum(OBmeankeep(:,1:numEvok),2)./(numEvok*Twin);
 %--- divide into diff stimuli ---
PC_spon1=sum(PCmnStim1(:,numEvok+1:numTot),2)./(numSpon*Twin);
PC_spon2=sum(PCmnStim2(:,numEvok+1:numTot),2)./(numSpon*Twin);
PC_od1=sum(PCmnStim1(:,1:numEvok),2)./(numEvok*Twin);
PC_od2=sum(PCmnStim2(:,1:numEvok),2)./(numEvok*Twin);
OB_spon1=sum(OBmnStim1(:,numEvok+1:numTot),2)./(numSpon*Twin);
OB_spon2=sum(OBmnStim2(:,numEvok+1:numTot),2)./(numSpon*Twin);
OB_od1=sum(OBmnStim1(:,1:numEvok),2)./(numEvok*Twin);
OB_od2=sum(OBmnStim2(:,1:numEvok),2)./(numEvok*Twin);

max_x=max([max([PC_sponBoth PC_odBoth]) max([OB_sponBoth OB_odBoth])]);
x=(0:.1: max_x)';
figure
hold on
plot(PC_sponBoth,PC_odBoth,'b*')
plot(OB_sponBoth,OB_odBoth,'ro')
plot(x,x,'k-.')
set(gca,'FontSize',18)
set(gca,'XLim',[0 max([max(OB_sponBoth) max(PC_sponBoth)])+1])
xlabel('Spontaneous Firing Rate (Hz)')
ylabel('Stim-Evoked Firing Rate (Hz)')
legend('PC','OB')


