function [ nspikes,counts, extracounts ] = extract_counts( spikes, UIDs, winT )
% extract_counts: turn spikes, UIDs into windowed spike counts
%
%   INPUTS:   spikes:     array of spikes
%             UIDs:       corresponding unit IDs
%             winT:       time window over which to count spikes
%
%   OUTPUTS:  counts:     array of counts
%             nspikes:    total spikes for each unit (can help to exclude
%                             non-spiking units)
%             extracounts:  list of spikes that occur after last winT
%             

%%% ASSUME minT = 0
maxT   = max(spikes);

maxID  = max(UIDs);


nwin  = round(maxT/winT);


counts   = zeros(nwin,maxID,'uint8');
nspikes  = zeros(maxID,1);

% Keep counts after last time window
extracounts = zeros(1,maxID);



for k=1:maxID
    rightID    = find(UIDs==k);
    
    rightSpk   = spikes(rightID);
    if(length(rightSpk)>1) %only count spikes separated by 0.1 ms
        rightSpk = rightSpk([1==1;diff(rightSpk)>.0001]); 
    end
    nspikes(k)=length(rightSpk);
    
   % [rightSpk double(UIDs(rightID))]
    
    if (~isempty(rightSpk))
        for k1=1:nwin
            %startT   = (k1-1)*winT; %never used b/c spikes are dynamically removed
            endT     = k1*winT;
            
            % First spike after end of interval
            lastind    = find(rightSpk > endT,1); 
            
            %
            if (isempty(lastind))
                % Any remaining spikes are in THIS interval
                % Count them and move on to 
                counts(k1,k) = length(rightSpk);
                break;
            elseif (lastind==1)
                % No spikes in this interval
                counts(k1,k) = 0;
            else
                % Any intervening spikes belong in this window
                counts(k1,k) = lastind-1;
                rightSpk = rightSpk(lastind:end);
            end
        end
        % Count any spikes AFTER last time window
        if (~isempty(rightSpk))
            extracounts(k) = length(find(rightSpk > nwin*winT));
        end
         
    end
        
end

