function [inds_before_drop,inds_before_rise,max_val] = findMaxCrosses(vals,threshold)

min_val = min(vals);
max_val = min_val+threshold;
changes = diff(vals>max_val);

inds_before_drop = find(changes==-1); % times the profile dipped below the threshold, i.e. where the min should par value should be
inds_before_rise = find(changes==1); % times the profile rose above the threshold, i.e. where the max should par value should be
