function [stim_trials,nEpoch,nPC,nOB]=processDataSelect(Twin,WinOvrLap,whichDatS)
%used in get_firing_rates.m and 
%relies on exclude_bad_units.m, extract_counts.m

if(whichDatS==0)
    infile = 'data040515_exp1+2.mat';
else
    infile = 'data040515_exp3+4.mat';
end

load(infile);

maxT  = ceil(max(max(OBspiketimes),max(aPCspiketimes)));



nPC  = max(aPCuid);
nOB  = max(OBuid);


%%%%% Firing rates
aPC_fr = zeros(nPC,1);
OB_fr  = zeros(nOB,1);

for k=1:nPC
    aPC_fr(k) = sum(aPCuid==k)/maxT;
end
for k=1:nOB
    OB_fr(k)  = sum(OBuid==k)/maxT;
end

avgfr_PC = sum(aPC_fr)/double(nPC);
avgfr_OB = sum(OB_fr)/double(nOB);

% Check:
%disp(['Average firing rate of PC: ' double(sum(aPC_fr)/double(nPC))])
%disp(['Average firing rate of OB: ' sum(OB_fr)/double(nOB)]) 



% Decide on thresholds for exclusion
valid_fr_LB  = 5./(maxT-stimtimes(1));
valid_fr_UB  = 49.;

% dont count anything before stimtimes(1); throwout rates\notin[valid_fr_LB,valid_fr_UB]
%i_aftFirstStimOB=OBspiketimes>stimtimes(1); %logical vector
[OBspiketimes, OBuid, OB_unit_map]=exclude_bad_units(OBspiketimes,OBuid,...
    valid_fr_LB,valid_fr_UB,maxT);
nOB = length(OB_unit_map);

%i_aftFirstStimPC=aPCspiketimes>stimtimes(1); %logical vector
[aPCspiketimes, aPCuid, PC_unit_map]=exclude_bad_units(aPCspiketimes,aPCuid,...
    valid_fr_LB,valid_fr_UB,maxT);
nPC = length(PC_unit_map);

%%%% Chop up data into trials
time_epochs  = [0 stimtimes maxT];

epoch_length = diff(time_epochs); 
nEpoch =length(epoch_length);

% Will include period before first stim
stim_trials = cell(nEpoch,1);
%PC_trials = cell(nEpoch,1);

%% Order the spikes
[OBsp_ord,I_OB]=sort(OBspiketimes,'ascend');
OBuid_ord = OBuid(I_OB);
[PCsp_ord,I_PC]=sort(aPCspiketimes,'ascend');
PCuid_ord = aPCuid(I_PC);


%% First set up a structure
for k=1:nEpoch
    stim_trials{k} = struct('OBspikes',[],'OBuid',[],...
                            'PCspikes',[],'PCuid',[],...
                            'startT',0,'endT',0,'stim_at_startT',0);
end

%% Enter some basic information
for k=1:nEpoch
    stim_trials{k}.startT = time_epochs(k);
    stim_trials{k}.endT   = time_epochs(k+1);
    if (k>1)
        stim_trials{k}.stim_at_startT = 1;
        % Already set to zero by default
    end
end
        

%% OB spike times
lastind = 0;
for k=1:nEpoch
   firstind = lastind+1; 
   if (k==nEpoch)
       lastind = length(OBsp_ord);
   else
       lastind  = find(OBsp_ord > time_epochs(k+1),1)-1;
   end
   [firstind lastind]
   
   % Record spikes relative to epoch start
   stim_trials{k}.OBspikes = OBsp_ord(firstind:lastind)-stim_trials{k}.startT;
   stim_trials{k}.OBuid    = OBuid_ord(firstind:lastind);
end
   
%% PC spike times
lastind = 0;
for k=1:nEpoch
   firstind = lastind+1; 
   if (k==nEpoch)
       lastind = length(PCsp_ord);
   else
       lastind  = find(PCsp_ord > time_epochs(k+1),1)-1;
   end
   [firstind lastind]
   stim_trials{k}.PCspikes = PCsp_ord(firstind:lastind)-stim_trials{k}.startT;
   stim_trials{k}.PCuid    = PCuid_ord(firstind:lastind);
end


% Spike counts
if(WinOvrLap==1)
    for k=1:nEpoch
        [nsp,cnts,Xcnts] = extractCountsOvrlap(stim_trials{k}.OBspikes,stim_trials{k}.OBuid,Twin);
        stim_trials{k}.OBcounts = cnts;
        
        [nsp,cnts,Xcnts] = extractCountsOvrlap(stim_trials{k}.PCspikes,stim_trials{k}.PCuid,Twin);
        stim_trials{k}.PCcounts = cnts;
        
    end
else %original; disjoint windows
    for k=1:nEpoch
        [nsp,cnts,Xcnts] = extract_counts(stim_trials{k}.OBspikes,stim_trials{k}.OBuid,Twin);
        % Check counts
        cksum = sum(abs(nsp-(sum(cnts)' + Xcnts')));
        if (~cksum)
            % Ok!
            stim_trials{k}.OBcounts = cnts;
        end
        
        [nsp,cnts,Xcnts] = extract_counts(stim_trials{k}.PCspikes,stim_trials{k}.PCuid,Twin);
        
        % Check counts
        cksum = sum(abs(nsp-(sum(cnts)' + Xcnts')));
        if (~cksum)
            % Ok!
            stim_trials{k}.PCcounts = cnts;
        end
    end
end



