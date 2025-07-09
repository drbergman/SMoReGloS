function profile = trimProfile(profile,threshold)

% This function takes in a full profile [npars+1,k] and a threshold and
% trims the profile to just the part below the threshold
% it uses linear interpolation to achieve this

%% get indices at which the profile is first below and then last below the threshold
[ind_first_below,ind_last_below] = getIndicesJustBelowThreshold(profile(end,:),threshold);
if numel(ind_first_below)>1 || numel(ind_last_below)>1
    error("This crosses below the threshold multiple times.")
end

max_val = min(profile(end,:)) + threshold;

if ind_first_below==1
    L = zeros(size(profile,1),0); % do not interpolate to the left because the first one is the first below the threshold
else
    slope = diff(profile(:,ind_first_below+(-1:0)),1,2);
    d = (max_val - profile(end,ind_first_below)) / diff(profile(end,ind_first_below+(-1:0)));
    L = profile(:,ind_first_below) + d * slope;
end

if ind_last_below==size(profile,2)
    R = zeros(size(profile,1),0); % do not interpolate to the right because the last one is the last below the threshold
else
    slope = diff(profile(:,ind_last_below+(0:1)),1,2);
    d = (max_val - profile(end,ind_last_below)) / diff(profile(end,ind_last_below+(0:1)));
    R = profile(:,ind_last_below) + d * slope;
end

profile = [L,profile(:,ind_first_below:ind_last_below),R];
