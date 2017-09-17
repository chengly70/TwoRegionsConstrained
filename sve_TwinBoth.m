% Script to save correl values as Twin varies; combining BOTH datasets
% using processDataSelect.m

Twin_vec=[0.05; 0.1; 0.2; 0.25; 0.4; 0.5; 1; 2];
Ltw=length(Twin_vec);

%variables to save; same naming convention as [mn/std][Cov/Cr][OB/PC/PB]_[spon/od][Both]
mnCovOB_sponBoth=zeros(Ltw,1);
stdCovOB_sponBoth=zeros(Ltw,1);
mnCovOB_odBoth=zeros(Ltw,1);
stdCovOB_odBoth=zeros(Ltw,1);
mnCovPC_sponBoth=zeros(Ltw,1);
stdCovPC_sponBoth=zeros(Ltw,1);
mnCovPC_odBoth=zeros(Ltw,1);
stdCovPC_odBoth=zeros(Ltw,1);
mnCovPB_sponBoth=zeros(Ltw,1);
stdCovPB_sponBoth=zeros(Ltw,1);
mnCovPB_odBoth=zeros(Ltw,1);
stdCovPB_odBoth=zeros(Ltw,1);
mnCrOB_sponBoth=zeros(Ltw,1);
stdCrOB_sponBoth=zeros(Ltw,1);
mnCrOB_odBoth=zeros(Ltw,1);
stdCrOB_odBoth=zeros(Ltw,1);
mnCrPC_sponBoth=zeros(Ltw,1);
stdCrPC_sponBoth=zeros(Ltw,1);
mnCrPC_odBoth=zeros(Ltw,1);
stdCrPC_odBoth=zeros(Ltw,1);
mnCrPB_sponBoth=zeros(Ltw,1);
stdCrPB_sponBoth=zeros(Ltw,1);
mnCrPB_odBoth=zeros(Ltw,1);
stdCrPB_odBoth=zeros(Ltw,1);
%-- adding in Fano Factor and Var ---
mnFFpc_sponBoth=zeros(Ltw,1);
stdFFpc_sponBoth=zeros(Ltw,1);
mnFFpc_odBoth=zeros(Ltw,1);
stdFFpc_odBoth=zeros(Ltw,1);
mnFFob_sponBoth=zeros(Ltw,1);
stdFFob_sponBoth=zeros(Ltw,1);
mnFFob_odBoth=zeros(Ltw,1);
stdFFob_odBoth=zeros(Ltw,1);
mnVrpc_sponBoth=zeros(Ltw,1);
stdVrpc_sponBoth=zeros(Ltw,1);
mnVrpc_odBoth=zeros(Ltw,1);
stdVrpc_odBoth=zeros(Ltw,1);
mnVrob_sponBoth=zeros(Ltw,1);
stdVrob_sponBoth=zeros(Ltw,1);
mnVrob_odBoth=zeros(Ltw,1);
stdVrob_odBoth=zeros(Ltw,1);

% -- main loop; loop over Twin_vec --
for twInd=1:Ltw
    VecCrPC_sponBoth=[]; %used to store data in all experiments
    VecCrPC_odBoth=[];
    VecCrOB_sponBoth=[];
    VecCrOB_odBoth=[];
    VecCrPB_sponBoth=[];
    VecCrPB_odBoth=[];
    
    VecCvPC_sponBoth=[];
    VecCvPC_odBoth=[];
    VecCvOB_sponBoth=[];
    VecCvOB_odBoth=[];
    VecCvPB_sponBoth=[];
    VecCvPB_odBoth=[];
    
    FFpc_odBoth=[];
    FFob_odBoth=[];
    FFpc_sponBoth=[];
    FFob_sponBoth=[];
    Vrpc_odBoth=[];
    Vrob_odBoth=[];
    Vrpc_sponBoth=[];
    Vrob_sponBoth=[];
    
    if(Twin_vec(twInd)~=2)
        WinOvrLap=1;
        numEvok = 2*(2/Twin_vec(twInd))-1; 
        numSpon = 2*(28/Twin_vec(twInd)); %first window can technically be in 2sec 'evoked' time
    else
        WinOvrLap=0;
        numEvok = (2/Twin_vec(twInd)); %assuming non-overlapping windows
        numSpon = (28/Twin_vec(twInd)); %discard data after 30 sec
    end
    
numTot=numEvok+numSpon;

for loopDataSet=1:2
% run processDataSelect.m
whichDatS=loopDataSet-1; %0 for Exp1+2, 1 for Exp3+4
[stim_trials,nEpoch,nPC,nOB]=processDataSelect(Twin_vec(twInd),WinOvrLap,whichDatS);

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

% Calc. Covariances (and correl) 
covPC_sponBoth=zeros(nPC,nPC); %a bit wasteful b/c symmetric
covPC_odBoth=zeros(nPC,nPC);
avgPC_sponBoth=zeros(nPC,1);
avgPC_odBoth=zeros(nPC,1);
varPC_sponBoth=zeros(nPC,1);
varPC_odBoth=zeros(nPC,1);
crrPC_sponBoth=zeros(nPC,nPC);
crrPC_odBoth=zeros(nPC,nPC);
covOB_sponBoth=zeros(nOB,nOB); %a bit wasteful b/c symmetric
covOB_odBoth=zeros(nOB,nOB);
avgOB_sponBoth=zeros(nOB,1);
avgOB_odBoth=zeros(nOB,1);
varOB_sponBoth=zeros(nOB,1);
varOB_odBoth=zeros(nOB,1);
crrOB_sponBoth=zeros(nOB,nOB);
crrOB_odBoth=zeros(nOB,nOB);
covPB_sponBoth=zeros(nPC,nOB); %a bit wasteful b/c symmetric
covPB_odBoth=zeros(nPC,nOB);
crrPB_sponBoth=zeros(nPC,nOB);
crrPB_odBoth=zeros(nPC,nOB);
%--odor(od), both --
%  PC-PC
tmpPCod=PCcount_allTrial_allCell(1:numEvok,:,i_PCs'); %get all data
tmpPCod=reshape(tmpPCod,numEvok*20,nPC); %matrix: numEvok*20 by nPC
avgPC_odBoth=sum(tmpPCod)'./(numEvok*20); %nPC by 1
tmpPCod=tmpPCod-repmat(avgPC_odBoth',numEvok*20,1); %centered
tmp2=tmpPCod; %for PC,PC cov/corr Evok
covPC_odBoth=(tmp2'*tmpPCod)./(numEvok*20-1); %unbiased estim. Cov
varPC_odBoth=diag(covPC_odBoth); %var spike counts, nPC by 1
crrPC_odBoth=covPC_odBoth./sqrt(varPC_odBoth*varPC_odBoth');
avg0s=avgPC_odBoth;
avg0s(avg0s==0)=Inf; %pad Fano with 0s
FFpc_odBoth=[FFpc_odBoth; varPC_odBoth./avg0s]; %forming FF here
Vrpc_odBoth=[Vrpc_odBoth; varPC_odBoth]; %forming Vr here
%  OB-OB
tmpOBod=OBcount_allTrial_allCell(1:numEvok,:,i_OBs'); %get all data
tmpOBod=reshape(tmpOBod,numEvok*20,nOB); %matrix: numEvok*20 by nOB
avgOB_odBoth=sum(tmpOBod)'./(numEvok*20); %nOB by 1
tmpOBod=tmpOBod-repmat(avgOB_odBoth',numEvok*20,1); %centered
tmp2=tmpOBod; %for OB,OB cov/corr Evok
covOB_odBoth=(tmp2'*tmpOBod)./(numEvok*20-1); %unbiased estim. Cov
varOB_odBoth=diag(covOB_odBoth); %var spike counts, nOB by 1
crrOB_odBoth=covOB_odBoth./sqrt(varOB_odBoth*varOB_odBoth');
avg0s=avgOB_odBoth;
avg0s(avg0s==0)=Inf; %pad Fano with 0s
FFob_odBoth=[FFob_odBoth; varOB_odBoth./avg0s]; %forming FF here
Vrob_odBoth=[Vrob_odBoth; varOB_odBoth]; %forming Var here

%  PC-OB
covPB_odBoth=(tmpPCod'*tmpOBod)./(numEvok*20-1);
crrPB_odBoth=covPB_odBoth./sqrt(varPC_odBoth*varOB_odBoth');
%--spont, both --
%  PC-PC
tmpPCsp=PCcount_allTrial_allCell(numEvok+1:numTot,:,i_PCs'); %get all data
tmpPCsp=reshape(tmpPCsp,numSpon*20,nPC); %matrix: numSpon*20 by nPC
avgPC_sponBoth=sum(tmpPCsp)'./(numSpon*20); %nPC by 1
tmpPCsp=tmpPCsp-repmat(avgPC_sponBoth',numSpon*20,1); %centered
tmp2=tmpPCsp; %for PC,PC cov/corr Evok
covPC_sponBoth=(tmp2'*tmpPCsp)./(numSpon*20-1); %unbiased estim. Cov
varPC_sponBoth=diag(covPC_sponBoth); %var spike counts, nPC by 1
crrPC_sponBoth=covPC_sponBoth./sqrt(varPC_sponBoth*varPC_sponBoth');
avg0s=avgPC_sponBoth;
avg0s(avg0s==0)=Inf; %pad Fano with 0s
FFpc_sponBoth=[FFpc_sponBoth; varPC_sponBoth./avg0s]; %store Fano here
Vrpc_sponBoth=[Vrpc_sponBoth; varPC_sponBoth]; %store Vr here

%  OB-OB
tmpOBsp=OBcount_allTrial_allCell(numEvok+1:numTot,:,i_OBs'); %get all data
tmpOBsp=reshape(tmpOBsp,numSpon*20,nOB); %matrix: numSpon*20 by nOB
avgOB_sponBoth=sum(tmpOBsp)'./(numSpon*20); %nOB by 1
tmpOBsp=tmpOBsp-repmat(avgOB_sponBoth',numSpon*20,1); %centered
tmp2=tmpOBsp; %for OB,OB cov/corr Evok
covOB_sponBoth=(tmp2'*tmpOBsp)./(numSpon*20-1); %unbiased estim. Cov
varOB_sponBoth=diag(covOB_sponBoth); %var spike counts, nOB by 1
crrOB_sponBoth=covOB_sponBoth./sqrt(varOB_sponBoth*varOB_sponBoth');
avg0s=avgOB_sponBoth;
avg0s(avg0s==0)=Inf; %pad Fano with 0s
FFob_sponBoth=[FFob_sponBoth; varOB_sponBoth./avg0s]; %store Fano here
Vrob_sponBoth=[Vrob_sponBoth; varOB_sponBoth]; %store Vr here
%  PC-OB
covPB_sponBoth=(tmpPCsp'*tmpOBsp)./(numSpon*20-1);
crrPB_sponBoth=covPB_sponBoth./sqrt(varPC_sponBoth*varOB_sponBoth');

% --- transform corr into vectors & plot results ---
idPC1=[];
idPC2=[];
for j=1:(nPC-1)
    idPC1=[idPC1; (1:j)'];
    idPC2=[idPC2; (j+1)*ones(j,1)];
end
ind_UpTri_PC=sub2ind([nPC nPC],idPC1,idPC2); %indices upper triang Corr/Cov matrix
idOB1=[];
idOB2=[];
for j=1:(nOB-1)
    idOB1=[idOB1; (1:j)'];
    idOB2=[idOB2; (j+1)*ones(j,1)];
end
ind_UpTri_OB=sub2ind([nOB nOB],idOB1,idOB2); %indices upper triang Corr/Cov matrix

%--- forming correl vectors to be saved ----
VecCrPC_sponBoth0=crrPC_sponBoth(ind_UpTri_PC);
VecCrPC_sponBoth=[VecCrPC_sponBoth; VecCrPC_sponBoth0(~isnan(VecCrPC_sponBoth0))]; %throw out 0 vars (NaN)
    VecCrPC_odBoth0=crrPC_odBoth(ind_UpTri_PC);
    VecCrPC_odBoth=[VecCrPC_odBoth; VecCrPC_odBoth0(~isnan(VecCrPC_odBoth0))]; %throw out 0 vars (NaN)
VecCrOB_sponBoth0=crrOB_sponBoth(ind_UpTri_OB);
VecCrOB_sponBoth=[VecCrOB_sponBoth; VecCrOB_sponBoth0(~isnan(VecCrOB_sponBoth0))]; %throw out 0 vars (NaN)
    VecCrOB_odBoth0=crrOB_odBoth(ind_UpTri_OB);
    VecCrOB_odBoth=[VecCrOB_odBoth; VecCrOB_odBoth0(~isnan(VecCrOB_odBoth0))]; %throw out 0 vars (NaN)
% PC-OB corre
VecCrPB_sponBoth0=reshape(crrPB_sponBoth,nPC*nOB,1);
VecCrPB_sponBoth=[VecCrPB_sponBoth; VecCrPB_sponBoth0(~isnan(VecCrPB_sponBoth0))]; %throw out 0 vars (NaN)
    VecCrPB_odBoth0=reshape(crrPB_odBoth,nPC*nOB,1);
    VecCrPB_odBoth=[VecCrPB_odBoth; VecCrPB_odBoth0(~isnan(VecCrPB_odBoth0))]; %throw out 0 vars (NaN)
%--- forming Cov vectors ----
VecCvPC_sponBoth=[VecCvPC_sponBoth; covPC_sponBoth(ind_UpTri_PC)];
    VecCvPC_odBoth=[VecCvPC_odBoth; covPC_odBoth(ind_UpTri_PC)];    
VecCvOB_sponBoth=[VecCvOB_sponBoth; covOB_sponBoth(ind_UpTri_OB)];
    VecCvOB_odBoth=[VecCvOB_odBoth; covOB_odBoth(ind_UpTri_OB)];
% PC-OB corre
VecCvPB_sponBoth=[VecCvPB_sponBoth; reshape(covPB_sponBoth,nPC*nOB,1)];
    VecCvPB_odBoth=[VecCvPB_odBoth; reshape(covPB_odBoth,nPC*nOB,1)];

end %end for loopDataSet
    
%save the results; combing both datasets
mnCrPC_sponBoth(twInd)=mean(VecCrPC_sponBoth); %store mean
stdCrPC_sponBoth(twInd)=std(VecCrPC_sponBoth); %store std
mnCrPC_odBoth(twInd)=mean(VecCrPC_odBoth); %store mean
stdCrPC_odBoth(twInd)=std(VecCrPC_odBoth); %store std
mnCrOB_sponBoth(twInd)=mean(VecCrOB_sponBoth); %store mean
stdCrOB_sponBoth(twInd)=std(VecCrOB_sponBoth); %store std
mnCrOB_odBoth(twInd)=mean(VecCrOB_odBoth); %store mean
stdCrOB_odBoth(twInd)=std(VecCrOB_odBoth); %store std
mnCrPB_sponBoth(twInd)=mean(VecCrPB_sponBoth); %store mean
stdCrPB_sponBoth(twInd)=std(VecCrPB_sponBoth); %store std
mnCrPB_odBoth(twInd)=mean(VecCrPB_odBoth); %store mean
stdCrPB_odBoth(twInd)=std(VecCrPB_odBoth); %store std
mnCovPC_sponBoth(twInd)=mean(VecCvPC_sponBoth); %store mean
stdCovPC_sponBoth(twInd)=std(VecCvPC_sponBoth); %store std
mnCovPC_odBoth(twInd)=mean(VecCvPC_odBoth); %store mean
stdCovPC_odBoth(twInd)=std(VecCvPC_odBoth); %store std
mnCovOB_sponBoth(twInd)=mean(VecCvOB_sponBoth); %store mean
stdCovOB_sponBoth(twInd)=std(VecCvOB_sponBoth); %store std
mnCovOB_odBoth(twInd)=mean(VecCvOB_odBoth); %store mean
stdCovOB_odBoth(twInd)=std(VecCvOB_odBoth); %store std
mnCovPB_sponBoth(twInd)=mean(VecCvPB_sponBoth); %store mean
stdCovPB_sponBoth(twInd)=std(VecCvPB_sponBoth); %store std
mnCovPB_odBoth(twInd)=mean(VecCvPB_odBoth); %store mean
stdCovPB_odBoth(twInd)=std(VecCvPB_odBoth); %store std
%-- fano factor --
mnFFpc_odBoth(twInd)=mean(FFpc_odBoth); %store mnFano here
stdFFpc_odBoth(twInd)=std(FFpc_odBoth); %store stdFano here
mnFFob_odBoth(twInd)=mean(FFob_odBoth); %store mnFano here
stdFFob_odBoth(twInd)=std(FFob_odBoth); %store stdFano here
mnFFpc_sponBoth(twInd)=mean(FFpc_sponBoth);
stdFFpc_sponBoth(twInd)=std(FFpc_sponBoth);
mnFFob_sponBoth(twInd)=mean(FFob_sponBoth);
stdFFob_sponBoth(twInd)=std(FFob_sponBoth);
% --- var spike cnts ---
mnVrpc_odBoth(twInd)=mean(Vrpc_odBoth); %store mnFano here
stdVrpc_odBoth(twInd)=std(Vrpc_odBoth); %store stdFano here
mnVrob_odBoth(twInd)=mean(Vrob_odBoth); %store mnFano here
stdVrob_odBoth(twInd)=std(Vrob_odBoth); %store stdFano here
mnVrpc_sponBoth(twInd)=mean(Vrpc_sponBoth);
stdVrpc_sponBoth(twInd)=std(Vrpc_sponBoth);
mnVrob_sponBoth(twInd)=mean(Vrob_sponBoth);
stdVrob_sponBoth(twInd)=std(Vrob_sponBoth);

end %end for twInd=1:Ltw


%%
save CoVar_Both mnC* stdC* Twin_vec mnFF* stdFF* mnVr* stdVr*

