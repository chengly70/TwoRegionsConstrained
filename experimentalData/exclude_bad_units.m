function [ spiketimes_out, UIDs_out, unit_map ] = exclude_bad_units( spiketimes, UIDs, ...
    valid_fr_LB, valid_fr_UB, maxT )
%  
%    exclude_bad_units:     remove units with firing rates below or above
%                           thresholds. Reorder unit numbers 
%    INPUTS:           spiketimes
%
    nUnits   = max(UIDs);
    fr_array = zeros(nUnits, 1);
    
    for k=1:nUnits
        fr_array(k) = sum(UIDs==k)/maxT;
    end
    
    Units_ok   = find(fr_array >= valid_fr_LB & fr_array <= valid_fr_UB);   
    Units_bad  = setdiff(1:nUnits,Units_ok);
    
    % Map from original -> new indices
    unit_map = [Units_ok [1:length(Units_ok)]'];

    for k=1:length(Units_bad)
        UIDs(find(UIDs==Units_bad(k)))=0;
    end

    for k=1:length(Units_ok)
        UIDs(find(UIDs==Units_ok(k)))=unit_map(k,2);
    end

    % Now remove "null" units
    spiketimes_out = spiketimes(find(UIDs));
    UIDs_out       = UIDs(find(UIDs));
    

end

