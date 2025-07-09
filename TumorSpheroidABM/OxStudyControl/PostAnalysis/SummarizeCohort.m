% This script will summarize my ABM output in a manner that will allow it
% to be used in the SMoRe-verse. This will serve as a first draft of a
% template for others to use to plug in their ABM to the SMoRe-verse.

% The output of this script is a new file called summary.mat with the
% following variables:
%   * t: vector of time values of data
%   * cohort_size: size of the cohort run. Expected use: grid sampling of ABM
%       parameters creates a [n1, n2, n3, ...] array of ABM parameter
%       vectors where each element corresponds to a node on the grid
%   * par_names: 1-D string vector of names of ABM parameters corresponding
%       to the dimensions in cohort_size
%   * vals: 1-D cell array of parameter values that correspond to each
%       dimension of the grid sample of ABM parameter space used for the 
%       cohort. vals{i} = [p_i_1, p_i_2, ..., p_i_ni] the ni values of ABM
%       parameter p_i that are varied along the first dimension of the
%       cohort
%   * nsamps_per_parameter_vector: number of ABM simulations run for each
%       cohort entry
%   * C: a 1-D cell array of containing data for the conditions under which ABM
%       simulations were run. Typical use case: running an ABM under
%       different dose concentrations of a therapeutic drug to match
%       experimentally-observed dose-response curves.
%   * n_conditions: number of conditions = numel(C) = length(C)
%   * n_time_series: number of features that are tracked over time in D below
%   * D: a struct array of the (D)ata. Each element corresponds to one
%       element of the cohort. For 1-D fields, the rows (first dimension)
%       corresponds to the time. For 2-D fields, the pages (third 
%       dimension) correspond to the time. The remaining dimensions
%       correspond to the feature. Of the following fields that the SMoRe-
%       verse recognizes, A must be included and either S or C must also be
%       included. D and Q are currently unused in the SMoRe-verse.
%       * A: average values of features, 1 feature per column
%       * S: standard deviation of features, 1 feature per column
%       * C: covariance matrix of features, D.C(:,:,i) is the covariance of
%          the features at time t(i)
%       * D: determinant of covariant matrix, D.D(i) is the determinant of
%          the D.C(:,:,i)
%       * Q: the inverse of the covariant matrix, D.Q(:,:,i) is the inverse
%           of D.C(:,:,i)

% a script to summarize cohort data. not sure why I originally created this
% as a function...but it seems more convenient to have this as a script.
% saves mean cell counts, phase counts, ode state variable counts, and all their
% SDs.


clearvars;
cohort_name = "cohort_2507081432";
% addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","nsamps_per_condition","cohort_size","lattice_parameters");
nsamps_per_parameter_vector = nsamps_per_condition;
n_conditions = 1;
C = {[]};
vals = {lattice_parameters.values};
par_names = strings(1,numel(lattice_parameters));
for i = 1:length(par_names)
    par_names(i) = lattice_parameters(i).path{end};
end
%%
for i = numel(ids):-1:1
    S = load(sprintf("../../data/sims/%s/output_final.mat",ids(i)));
    if size(S.tracked.phases,2)>4
        S.tracked.phases = reshape(S.tracked.phases,[],4,2);
        S.tracked.phases = sum(S.tracked.phases,3);
        % error("Not sure how we will count the arrested compartment in this.")
    end
    phase_count(:,:,:,i) = reshape(S.tracked.phases,[],2,2);
end

t_abm = S.tracked.t;
t_abm = round(1440*t_abm)/1440; % to make sure that the last time point is actually 3 days (not 3-eps() days)
nt_abm = length(t_abm);

t = [0;10;24;36;48;72]/24;
nt = length(t);


%% reshape and summarize
phase_count = reshape(phase_count,[nt_abm,2,2,size(ids)]);
phase_count_sampled = interp1(t_abm,phase_count,t);
state_vars = sum(phase_count_sampled,2);
state_vars = reshape(state_vars,nt,2,[],nsamps_per_parameter_vector);

if nsamps_per_parameter_vector==1
    temp_avg = state_vars;
    temp_std = zeros(size(state_vars));
else
    temp_avg = mean(state_vars,ndims(state_vars));
    temp_std = std(state_vars,[],ndims(state_vars));
end

temp_avg = reshape(temp_avg,nt,2,[]);
temp_std = reshape(temp_std,nt,2,[]);

%% create data struct and replace 0 SD with nonzero SD
% if isempty(gcp('nocreate'))
%     parpool("Threads")
% end
r_threshold = 1e-14;
for i = size(temp_avg,3):-1:1
    D(i).A = temp_avg(:,:,i);
    covs = zeros(2,2,nt);
    r = zeros(nt,1);
    % parfor ti = 1:nt
    for ti = 1:nt
        covs(:,:,ti) = cov(squeeze(state_vars(ti,:,i,:))');
        r(ti) = rcond(squeeze(covs(:,:,ti)));
    end

    bad_r = r < r_threshold;
    I = find(bad_r);
    for ii = length(I):-1:1
        if I(ii)==1 && ~bad_r(2) % if it is the first time point and the second is nonzero, use that
            covs(:,:,1) = covs(:,:,2);
            r(1) = r(2);
            bad_r(1) = bad_r(2);
        elseif I(ii)==nt && ~bad_r(nt-1) % similarly for the last time point
            covs(:,:,nt) = covs(:,:,nt-1);
            r(nt) = r(nt-1);
            bad_r(nt) = bad_r(nt-1);
        elseif all(I(ii)~=[1,nt]) % if it's in the middle...
            if ~any(bad_r(I(ii)+[-1,1])) % take the average of the SDs on either side
                covs(:,:,I(ii)) = 0.5*(covs(:,:,I(ii)-1)+covs(:,:,I(ii)+1));
                r(I(ii)) = rcond(squeeze(covs(:,:,I(ii))));
                bad_r(I(ii)) = r(I(ii)) < r_threshold;
                if bad_r(I(ii))
                    covs(:,:,I(ii)) = eye(2);
                    r(I(ii)) = 1;
                    bad_r(I(ii)) = false;
                end
            elseif ~bad_r(I(ii)-1) % if the one prior is nonzero, use that
                covs(:,:,I(ii)) = covs(:,:,I(ii)-1);
                r(I(ii)) = r(I(ii)-1);
                bad_r(I(ii)) = bad_r(I(ii)-1);
            elseif ~bad_r(I(ii)+1) % if the one after is nonzero, use that
                covs(:,:,I(ii)) = covs(:,:,I(ii)+1);
                r(I(ii)) = r(I(ii)+1);
                bad_r(I(ii)) = bad_r(I(ii)+1);
            else % if all else fails, set it to Identity to use the actual error
                covs(:,:,I(ii)) = eye(2);
                r(I(ii)) = 1;
                bad_r(I(ii)) = false;
            end
        else % then it is at the beginning or end with a zero next to it
            covs(:,:,I(ii)) = eye(2);
            r(I(ii)) = 1;
            bad_r(I(ii)) = false;
        end

        if bad_r(I(ii))
            covs(:,:,I(ii)) = eye(2);
            r(I(ii)) = 1;
            bad_r(I(ii)) = false;
        end

    end

    dets = zeros(nt,1);
    invs = zeros(2,2,nt);
    % parfor ti = 1:nt
    for ti = 1:nt
        dets(ti) = det(covs(:,:,ti));
        invs(:,:,ti) = inv(covs(:,:,ti));
    end
    D(i).C = covs;
    D(i).D = dets;
    D(i).Q = invs;
end

for i = size(temp_avg,3):-1:1
    [I,J] = find(temp_std(:,:,i)==0);
    for ii = 1:length(I)
        if I(ii)==1 && temp_std(2,J(ii),i)~=0 % if it is the first time point and the second is nonzero, use that
            temp_std(1,J(ii),i) = temp_std(2,J(ii),i);
        elseif I(ii)==nt && temp_std(nt-1,J(ii),i)~=0 % similarly for the last time point
            temp_std(nt,J(ii),i) = temp_std(nt-1,J(ii),i);
        elseif all(I(ii)~=[1,nt]) % if it's in the middle...
            if all(temp_std(I(ii)+[-1,1],J(ii),i)~=0) % take the average of the SDs on either side
                temp_std(I(ii),J(ii),i) = mean(temp_std(I(ii)+[-1,1],J(ii),i));
            elseif any(temp_std(I(ii)+[-1,1],J(ii),i)~=0) % if one is nonzero, use that
                temp_std(I(ii),J(ii),i) = max(temp_std(I(ii)+[-1,1],J(ii),i));
            else % if all else fails, set it to 1 to use the acutal error
                temp_std(I(ii),J(ii),i) = 1; 
            end
        else % then it is at the beginning or end with a zero next to it
            temp_std(I(ii),J(ii),i) = 1;
        end
    end
    D(i).S = temp_std(:,:,i);
end
D = reshape(D,[1,cohort_size]);


n_time_series = size(D(1).A,2);
save(sprintf("../../data/%s/summary.mat",cohort_name),"D","t","C","cohort_size","nsamps_per_parameter_vector","n_conditions","vals","n_time_series","par_names","-v7.3")