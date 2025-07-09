function [par_min,par_max] = getProfileBounds(profile,threshold)

% This will likely be phased out.

% takes in profile and returns the bounds on the parameter that keeps the
% profile within the threshold
% profile: a 2xn array with row 1 being parameter values and row 2 being
% goodness-of-fit values
% threshold: chi2inv value that determines difference between max permitted
% value and the min of the second row in profile
% this throws an error if the threshold is crossed too many times, which
% occurs when the optimization found different local minima or possibly if
% the profile has multiple local minima below the threshold

assert(size(profile,1)==2)

[inds_before_drop,inds_before_rise,max_val] = findMaxCrosses(profile(2,:),threshold);

if numel(inds_before_drop)>1 || numel(inds_before_rise)>1 % need to explore why so many crosses were found, unclear how to programmatically proceed
    error("Too many crosses found.")
end

if ~isempty(inds_before_rise) && ~isempty(inds_before_drop) && inds_before_drop>inds_before_rise
    error("Inverted profile shape that rises above the threshold then falls back below. Need to explore more.")
end

if isempty(inds_before_drop)
    par_min = profile(1,1);
else % numel(inds_before_drop)==1 % use linear interpolation to estimate where in this change the threshold is hit
    x = profile(1,inds_before_drop+(0:1));
    y = profile(2,inds_before_drop+(0:1));
    par_min = x(1) + diff(x)*(max_val-y(1))/diff(y);
end

if isempty(inds_before_rise)
    par_max = profile(1,end);
else % numel(inds_before_rise)==1 % use linear interpolation to estimate where in this change the threshold is hit
    x = profile(1,inds_before_rise+(0:1));
    y = profile(2,inds_before_rise+(0:1));
    par_max = x(1) + diff(x)*(max_val-y(1))/diff(y);
end
