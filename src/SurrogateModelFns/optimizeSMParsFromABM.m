function [P,fstar] = optimizeSMParsFromABM(files,sm,p0,lb,ub,optim_opts,weights,opts,raw_error_opts)

% THIS IS A USER-FACING FUNCTION

% Function to optimize SM parameters at various points in ABM parameter
% space.
%   * files: structure with fields that are strings with paths to certain
%       files:
%       * data: summarized data of the ABM simulations (see
%           SMoReParS/TumorSpheroidABM/OxStudyControl/PostAnalysis/SummarizeCohort.m
%           for explanation of the loaded variables
%       * previous_optim_file: checkpoint of this function
%   * p0: default initial parameter values for SM
%   * sm: structure of SM information that has the following fields:
%       * fn: function to simulate the SM
%       * opts: options to run with the SM. any options that are needed by the function; use this for
%           options that do NOT vary within one call to the objective function, e.g.
%           options that specify the version of the SM
%   * lb: lower bound of SM parameter values for optimization
%   * ub: upper bound of SM parameter values for optimization
%   * optim_opts: optimization options for fmincon
%   * weights: weights to be used across conditions, does not need to be
%       normalized
%   * opts: options for this function. See arguments block for the 
%       available options. supply as name-value pairs, e.g.
%       optimizeSMParsFromABM(...,force_serial=true)

arguments
    files struct
    sm struct
    p0 (:,1) double
    lb (:,1) double
    ub (:,1) double
    optim_opts struct
    weights (:,1) double
    opts.force_serial logical = false
    opts.n_starts {mustBeInteger} = 1
    opts.save_every_iter {mustBeInteger} = 1e300
    opts.save_every_sec double = Inf
    opts.checkpoint_filename string = sprintf("data/temp_profile_%02d",next_version_number("data/temp_profile_%02d"))

    % rawError opts
    % WARNING: These should be updated in tandem with the logical flow of
    %           rawError to ensure consistency between the two.
    raw_error_opts.assume_independent_time_series = true; % assume that the time series produced by the SM are independent (crazy, right? but it's what I had initially assumed, so this is the default value)
    raw_error_opts.only_use_z_scores = true; % whether to use the constant and SD terms from the LL for normal distributions, or (if false) just use the sum of z-scores
    raw_error_opts.report_as_error = true; % whether to report the value as an error for optimization purposes
    raw_error_opts.resample_t = []; % time points to resample at and compare with data, leave empty to not resample but to use time points from D
end

load(files.data,"t","D","C");

cohort_size = size(D,2:ndims(D)); % size of D is [n conditions, cohort_size]
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used
m = size(D,1); % number of conditions used

assert(isempty(C) || (numel(C)==m && numel(weights)==m)) % make sure that these match if there are any conditions

npars = numel(p0);
sz = [opts.n_starts,n_abm_vecs];
n_total = prod(sz);
if isfield(files,"previous_optim_file")
    load(files.previous_optim_file,"p_all","f_all");
    if size(f_all,1)>opts.n_starts % now asking for fewer initializations for optimization
        f_all = f_all(1:opts.n_starts,:);
        p_all = p_all(:,1:opts.n_starts,:);
    elseif size(f_all,1)<opts.n_starts % now asking for more initializations for optimization
        n_new = opts.n_starts-size(f_all,1);
        f_all = cat(1,f_all,zeros(n_new,n_abm_vecs));
        p_all = cat(2,p_all,zeros(npars,n_new,n_abm_vecs));
    end
    p_all = reshape(p_all,[npars,sz]);
    f_all = reshape(f_all,sz);
    ind_to_run = find(f_all(:)==0);
else
    p_all = zeros([npars,sz]);
    f_all = zeros(sz);
    ind_to_run = 1:n_total;
end
num_to_run = numel(ind_to_run);

if isfield(sm,"custom_raw_error_fn")
    rE_fn = @(p,d,c) sm.custom_raw_error_fn(sm,p,t,d,c,raw_error_opts);
else
    clear rawError
    rE_fn = @(p,d,c) rawError(sm,p,t,d,c,raw_error_opts);
end

last_save_time = tic;

data_dir = split(opts.checkpoint_filename,"/");
data_dir = join(data_dir(1:end-1),"/");
warning("off",'MATLAB:MKDIR:DirectoryExists')
mkdir(data_dir)
warning("on",'MATLAB:MKDIR:DirectoryExists')

if ~opts.force_serial
    if isempty(gcp('nocreate'))
        ppool = parpool("Processes");
    elseif isa(gcp,"parallel.ThreadPool")
        delete(gcp('nocreate'))
        ppool = parpool("Processes");
    else
        ppool = gcp;
    end

    %% This submits them all using parfeval
    FF(1:num_to_run) = parallel.FevalFuture;
    for i = 1:num_to_run
        [sample_ind,abm_par_ind] = ind2sub(sz,ind_to_run(i));
        if m==1
            d = D(abm_par_ind);
            F = @(p) rE_fn(p,d,C{1});
        else
            d = D(:,abm_par_ind);
            F = @(p) arrayfun(@(j) rE_fn(p,d(j),C{j}),1:m)*weights;
        end
        if sample_ind>1
            p0_si = min(ub,max(lb,p0.*exp(.5*randn(npars,1))));
        else
            p0_si = p0;
        end
        FF(i) = parfeval(ppool,@() fmincon(F,p0_si,[],[],[],[],lb,ub,[],optim_opts),2);
    end

    for i = 1:num_to_run
        [temp_idx,p_temp,f_temp] = fetchNext(FF);
        [sample_ind,abm_par_ind] = ind2sub(sz,ind_to_run(temp_idx));
        p_all(:,sample_ind,abm_par_ind) = p_temp;
        f_all(sample_ind,abm_par_ind) = f_temp;

        if mod(i,round(.01*num_to_run))==0
            fprintf("Finished %d of %d.\n",i,num_to_run)
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.checkpoint_filename,"p_all","f_all")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end

else
    p_all = reshape(p_all,[npars,sz]);
    for i = 1:num_to_run
        [sample_ind,abm_par_ind] = ind2sub(sz,ind_to_run(i));

        if sample_ind>1
            p0_si = min(ub,max(lb,p0.*exp(.5*randn(npars,1))));
        else
            p0_si = p0;
        end
        if m==1
            [p_all(:,sample_ind,abm_par_ind),f_all(sample_ind,abm_par_ind)] = fmincon(@(p) rE_fn(p,D(abm_par_ind),C{1}),p0_si,[],[],[],[],lb,ub,[],optim_opts);
        else
            [p_all(:,sample_ind,abm_par_ind),f_all(sample_ind,abm_par_ind)] = fmincon(@(p) arrayfun(@(j) rE_fn(p,D(j,abm_par_ind),C{j}),1:m)*weights,p0_si,[],[],[],[],lb,ub,[],optim_opts);
        end

        if mod(i,round(.01*num_to_run))==0
            fprintf("Finished %d of %d.\n",i,num_to_run)
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.checkpoint_filename,"p_all","f_all")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end


end

[fstar,best_case_inds] = min(f_all,[],1);
fstar = reshape(fstar,cohort_size);
P = zeros(npars,n_abm_vecs);
 
for i = 1:n_abm_vecs
    P(:,i) = p_all(:,best_case_inds(i),i);
end

P = reshape(P,[npars,cohort_size]);

end
