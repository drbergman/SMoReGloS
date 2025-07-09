function [x,y_mean,patch_coords] = patchPlotCoords(X,Y,input_opts)

% this will plot the mean and +/- 1 SD patch

% X is the vector/array of independent variable values. If an array, then
% the size(X,1) is the number of time points and size(X,2) is the number of
% samples

% Y is the array of dependent variables from n samples where size(Y,2)==n.

opts = defaultPatchPlotCoordsOptions;
if nargin==3 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

if isvector(X)
    X = reshape(X,[],1);
end

posdiff_or_nan = @(x) x>0 | isnan(x);
assert(size(X,1)==size(Y,1) || length(X)==size(Y,1)) % make sure there are equal observations of indepdent and dependent variables
for ci = 1:size(X,2)
    if isnan(X(end,ci)) % then make sure that all NaNs lead up to this
        assert(all(diff(find(isnan(X(:,ci))))==1))
    end
end
assert(all(posdiff_or_nan(diff(X,1,1)),'all')) % make sure all columns of X are increasing or nan
% assert(numel(unique(X(1,:)))==1) % make sure all columns of X start at same time point
% assert(numel(unique(X(end,:)))==1) % make sure all columns of X end at same time point

if size(X,2) > 1 % then each sample has own time series, normalize this
    x = unique(X(:));
    y = NaN(length(x),size(X,2));
    for ci = 1:size(X,2)
        ntemp = find(~isnan(X(:,ci)),1,"last");
        y(:,ci) = interp1(X(1:ntemp,ci),Y(1:ntemp,ci),x,"linear");
    end
else
    x = X;
    y = Y;
end

if opts.omit_nan==true
    y_mean = mean(y,2,"omitnan");
    nan_log = isnan(y_mean); % where all the values were NaN
    x(nan_log) = [];
    y(nan_log,:) = [];
    y_mean(nan_log) = [];
    if opts.split_sd
        y_nan_log = isnan(y);
        y_lower_log = ~y_nan_log & (y < y_mean);
        y_upper_log = ~y_nan_log & (y > y_mean);
        y_at_mean_log = ~y_nan_log & ~y_lower_log & ~y_upper_log;
        [ri_lower,~] = find(y_lower_log);
        [ri_upper,~] = find(y_upper_log);
        dy2 = (y - y_mean).^2;
        n_at_mean = sum(y_at_mean_log,2);
        n_lower = sum(y_lower_log,2)+0.5*n_at_mean;
        n_upper = sum(y_upper_log,2)+0.5*n_at_mean;
        sd_lower = sqrt(accumarray(ri_lower,dy2(y_lower_log),[size(y,1),1])./n_lower);
        sd_upper = sqrt(accumarray(ri_upper,dy2(y_upper_log),[size(y,1),1])./n_upper);
        y_sd = [sd_lower,sd_upper];
    else
        y_sd = repmat(std(y,[],2,"omitnan"),1,2);
    end
else
    y_mean = mean(y,2);
    if opts.split_sd
        y_lower_log = (y < y_mean);
        y_upper_log = (y > y_mean);
        y_at_mean_log = ~y_lower_log & ~y_upper_log;
        [ri_lower,~] = find(y_lower_log);
        [ri_upper,~] = find(y_upper_log);
        dy2 = (y - y_mean).^2;
        n_at_mean = sum(y_at_mean_log,2);
        n_lower = sum(y_lower_log,2)+0.5*n_at_mean;
        n_upper = sum(y_lower_log,2)+0.5*n_at_mean;
        sd_lower = sqrt(accumarray(ri_lower,dy2(y_lower_log),[size(y,1),1])/n_lower);
        sd_upper = sqrt(accumarray(ri_upper,dy2(y_upper_log),[size(y,1),1])/n_upper);
        y_sd = [sd_lower,sd_upper];
    else
        y_sd = repmat(std(y,[],2),1,2);
    end
end

patch_coords{1} = [x;flip(x)];
patch_coords{2} = [y_mean-y_sd(:,1);flip(y_mean+y_sd(:,2))];

end

function default_options = defaultPatchPlotCoordsOptions

default_options.omit_nan = false;
default_options.split_sd = false; % if false, use the normal standard deviation, otherwise split it into upper and lower SD

end