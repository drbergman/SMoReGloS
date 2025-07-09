function [I1,I2] = getProfileBoundaryInds(vals,threshold)

% This will likely be phased out.

% takes in profile vals and returns the bounds [I1,I2] on the indices that keeps 
% the profile within the threshold
% profile: a 1xn array of goodness-of-fit values
% threshold: chi2inv value that determines difference between max permitted
% value and the min of the second row in profile

[inds_before_drop,inds_before_rise,~] = findMaxCrosses(vals,threshold);

if numel(inds_before_drop)>1 || numel(inds_before_rise)>1 % need to explore why so many crosses were found, unclear how to programmatically proceed
    error("Too many crosses found.")
end

if ~isempty(inds_before_rise) && ~isempty(inds_before_drop) && inds_before_drop>inds_before_rise
    error("Inverted profile shape that rises above the threshold then falls back below. Need to explore more.")
end

if isempty(inds_before_drop) % then it started below the threshold and never temporarily increased beyond the threshold
    I1 = 1;
else % numel(inds_before_drop)==1 % the next index is the first to drop below the max value
    I1 = inds_before_drop + 1;
end

if isempty(inds_before_rise) % then it never rose above the threshold, can go all the way to the end
    I2 = length(vals);
else % numel(inds_before_rise)==1 % then this is actually the index to use
    I2 = inds_before_rise;
end
