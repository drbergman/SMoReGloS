function [inds_first_below,inds_last_below] = getIndicesJustBelowThreshold(vals,threshold)

% finds all indices that are the first (in a sequence from left-to-right)
% that are below the threshold
% similarly, finds the last that are below the threshold
[inds_before_drop,inds_last_below,max_val] = findMaxCrosses(vals,threshold);
inds_first_below = inds_before_drop+1; % this is then the index is first after dropping below the threshold
if vals(1)<=max_val % if the profile begins below the max_val, add this index to the list (note, 1 cannot occur above due to the +1 to indices >=1)
    inds_first_below = [1,inds_first_below];
end

if vals(end)<=max_val % if the profile ends below the max_val, add the final index to this list (note, n cannot already be in the list because diff returns a vector of length n-1)
    inds_last_below = [inds_last_below,length(vals)];
end
