function raw_error = getRawError(sm,p,tt,D,C,opts)

% THIS IS A USER-FACING FUNCTION

% a user-friendly wrapper for rawError that allows the user to supply
% name-value arugments, e.g. getRawError(..., report_as_error=false), to
% return the result of rawError

arguments
    sm
    p (:,1) double
    tt (:,1) double
    D
    C

    % WARNING: These should be updated in tandem with the logical flow of
    %           rawError to ensure consistency between the two.
    opts.assume_independent_time_series = true; % assume that the time series produced by the SM are independent (crazy, right? but it's what I had initially assumed, so this is the default value)
    opts.only_use_z_scores = true; % whether to use the constant and SD terms from the LL for normal distributions, or (if false) just use the sum of z-scores
    opts.report_as_error = true; % whether to report the value as an error for optimization purposes
    opts.resample_t = []; % time points to resample at and compare with data, leave empty to not resample but to use time points from D
end

raw_error = rawError(sm,p,tt,D,C,opts);