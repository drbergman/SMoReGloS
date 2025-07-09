function out = profileLikelihood_one_ind(profile_ind,pbest,tt,data,data_std,para_ranges,lb,ub,opts,threshold,save_all_pars)

% profiles each parameter in pbest. pbest is the best fit at the point. the
% time series to compare against is given in tt, data, and data_std.
% para_ranges sets bounds for profiling. lb and ub set bounds for
% optimization while another parameter is being profiled. opts carries
% fmincon options. threshold sets a cutoff for when profiling can be
% stopped. save_all_pars can be passed in as true to save all parameter
% values, not just the profiled parameter value

% In each direction of pbest(i) (left and right), the ith parameter changes
% at 1% of pbest(i) for 10 steps. If it is still within para_ranges(i,:)
% and the value of the output is within F(pbest)+threshold, then continue
% in that direction until one of those conditions fails.

addpath("../ODEFitting/")

if nargin < 10
    save_all_pars = false;
end
npars = numel(pbest);

if size(data,2)==2 % if two time series are given, then compare against both
    F = @(p) sum(((computeTimeSeries(p,tt) - data)./data_std).^2,'all'); % if the data has two phases
else % otherwise, compare total (summed) counts
    F = @(p) sum(((sum(computeTimeSeries(p,tt),2) - data)./data_std).^2,'all'); % if the data has just total cell count
end

for i = profile_ind

    Aeq = zeros(1,npars);
    Aeq(i) = 1;
    min_val = F(pbest);
    dxi = pbest(i)*0.01; % move at 1% of the best fit val
    for dir = [-1,1] % move left and right along this parameter dimension
        x0 = pbest;
        temp_val = zeros(1,10);
        if save_all_pars
            temp_par = zeros(npars,10);
        else
            temp_par = zeros(1,10);
        end
        for j = 1:10
            x0(i) = x0(i) + dir*dxi;
            [x0,temp_val(j)] = fmincon(F,x0,[],[],Aeq,x0(i),lb,ub,[],opts);
            if save_all_pars
                temp_par(:,j) = x0;
            else
                temp_par(j) = x0(i);
            end
        
        end
        while (x0(i) + dir*dxi >= para_ranges(i,1)) && (x0(i) + dir*dxi <= para_ranges(i,2)) && (temp_val(end) < min_val + threshold)
            x0(i) = x0(i) + dir*dxi;
            if save_all_pars
                temp_par(:,end+1) = x0;
            else
                temp_par(end+1) = x0(i);
            end
            [x0,temp_val(end+1)] = fmincon(F,x0,[],[],Aeq,x0(i),lb,ub,[],opts);
        end
        if dir==-1
            if save_all_pars
                out = [fliplr(temp_par),pbest';flip(temp_val),min_val];
            else
                out = [fliplr(temp_par),pbest(i);flip(temp_val),min_val];
            end
        else
            out = [out,[temp_par;temp_val]];
        end
    end

end

rmpath("../ODEFitting/")
