function profiles = performProfile(files,sm,profile_params,opts,profile_likelihood_opts,raw_error_opts)

% THIS IS A USER-FACING FUNCTION

% profiles each SM parameter for each "column" of the data array, Data
% array being reshaped to a 2D array.
%
% files: struct with fields
%   optimal_parameters: best fit SM parameters for each column of Data
%   data: the data in size [nconditions,cohort_size]
%   previous_profile_file (optional): partially completed profile
% objfn_constants: struct with fields that do vary by condition
%   fn: function used by rawError to generate SM output to compare against data
%   fn_opts: options to be used by fn for all conditions
%   weights: weights to use across conditions for objective function
% profile_params: struct with fields
%   lb: lower bound of SM parameters in calls to fmincon
%   ub: upper bound of SM parameters in calls to fmincon
%   A (optional): array to give A*x<=b constraint in fmincon
%   b (optional): vector to give A*x<=b constraint in fmincon
%   opts: options to pass in to fmincon
%   step_growth_factor (optional): factor to increase step size in Stage 2 of profiling
%   para_ranges: n_sm_pars x 2 array of extrema to be used in profiling for each parameter, i.e. profiles never force this parameter beyond this value
%   shrinking_factor: factor by which to shrink step sizes when approaching bounds in para_ranges
%   initial_step_prop: proportion of best parameter to move in Stage 1
%   smallest_par_step: minimum allowable step size for each parameter; Stage 2 halts when the step size is below this and the threshold is exceeded
%   min_num_steps: number of steps for each SM to take in Stage 1
%   threshold: chi2inv value for profile
%   secondary_step_factor: factor to increase an SM parameter step size when beginning Stage 2
%   step_growth_factor: factor to increase step size in Stage 2 after successfully extending profile
% opts (optional): struct with (all optional) fields
%   force_serial=true: logical to force serial computation of profiles
%   save_all_pars=true: logical to control whether to save all parameter values (true) or only the current profile parameter and goodness-of-fit value (false)
%   save_every_iter=Inf: how often (based on iterations) to save profile throughout (protects against errors and workers crashing)
%   save_every_sec=Inf: how often (based on seconds passed) to save profile throughout (protects against errors and workers crashing)
%   temp_profile_name: how often to save profile throughout (protects against errors and workers crashing)
%   profile_likelihood_options: a structure of options to be passed in to profileLikelihood.m with (all optional) fields
%       save_all_pars=true: logical to control whether to save all parameter values (true) or only the current profile parameter and goodness-of-fit value (false)
%       raw_error_opts=[]: if not empty, then a struct of options to be used in rawError.m

arguments
    files struct
    sm struct
    profile_params struct

    opts.force_serial logical = true
    opts.save_every_iter {mustBeInteger} = 1e300; % how often (based on iterations) to save profile throughout (protects against errors and workers crashing)
    opts.save_every_sec double = Inf; % how often (based on seconds passed) to save profile throughout (protects against errors and workers crashing)
    opts.temp_profile_name string = sprintf("data/temp_profile_%02d",next_version_number("data/temp_profile_%02d"));

    % profileLikelihood opts
    % WARNING: These should be updated in tandem with the logical flow of
    %           profileLikelihood to ensure consistency between the two.
    profile_likelihood_opts.save_all_pars logical = true;

    % rawError opts
    % WARNING: These should be updated in tandem with the logical flow of
    %           rawError to ensure consistency between the two.
    raw_error_opts.assume_independent_time_series logical = true; % assume that the time series produced by the SM are independent (crazy, right? but it's what I had initially assumed, so this is the default value)
    raw_error_opts.only_use_z_scores logical = true; % whether to use the constant and SD terms from the LL for normal distributions, or (if false) just use the sum of z-scores
    raw_error_opts.report_as_error logical = true; % whether to report the value as an error for optimization purposes
    raw_error_opts.resample_t double = []; % time points to resample at and compare with data, leave empty to not resample but to use time points from D
end

if ~isfield(profile_params,"A")
    profile_params.A = [];
    profile_params.b = [];
end

if isfield(files,"optimal_parameters")
    load(files.optimal_parameters,"P")
    if isequal(sm,struct())
        load(files.optimal_parameters,"sm")
    end
    if ~isfield(files,"data")
        temp = load(files.optimal_parameters,"files");
        if isfield(temp,"files") && isfield(temp.files,"data")
            load(temp.files.data,"t","D","C","cohort_size");
        end
    end
    if ~isfield(profile_params,"lb")
        temp = load(files.optimal_parameters,"lb","ub");
        profile_params.lb = temp.lb;
        profile_params.ub = temp.ub;
    end
    if ~isfield(profile_params,"opts")
        temp = load(files.optimal_parameters,"optim_opts");
        profile_params.opts = temp.optim_opts;
    end
    if ~isfield(profile_params,"weights")
        temp = load(files.optimal_parameters,"weights");
        profile_params.weights = temp.weights;
    end
else
    error("Rename this field in files from par_file --> optimal_parameters")
end

if isfield(files,"data")
    load(files.data,"t","D","C","cohort_size");
elseif isfield(files,"data_file")
    warning("Rename this field in files from data_file --> data")
    load(files.data_file,"t","D","C","cohort_size");
end
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used

P = reshape(P,size(P,1),[]);
npars = size(P,1);
t_start = tic;
if isfield(files,"previous_profile_file")
    load(files.previous_profile_file,"profiles")
    profiles = reshape(profiles,npars,n_abm_vecs);
    ind_to_run = find(any(cellfun(@isempty,profiles),1));
else
    profiles = cell(npars,n_abm_vecs);
    ind_to_run = 1:n_abm_vecs;
end

num_to_run = length(ind_to_run);
fprintf("Going to run (# SM Parameters) x (# ABM Parameters) = %d x %d = %d profiles.\n",npars,num_to_run,npars*num_to_run);

if ~isfield(profile_params,"step_growth_factor")
    profile_params.step_growth_factor = ones(npars,1);
end

last_save_time = tic;

% pL_fn = @(p,d) profileLikelihood(p,t,d,C,sm,profile_params,profile_likelihood_opts,raw_error_opts);

if ~opts.force_serial
    %% run in parallel
    FF(1:num_to_run) = parallel.FevalFuture;

    if isempty(gcp('nocreate'))
        ppool = parpool("Processes");
    else
        ppool = gcp;
    end
    for i = 1:num_to_run
        current_ind = ind_to_run(i);
        d = D(:,current_ind);
        p = P(:,current_ind);
        FF(i) = parfeval(ppool,@() profileLikelihood(p,t,d,C,sm,profile_params,profile_likelihood_opts,raw_error_opts),1);
    end
    fprintf("FevalQueue finished.\n")
    %% fetch profiles performed in parallel
    for i = 1:num_to_run

        [idx,temp] = fetchNext(FF);
        profiles(:,ind_to_run(idx)) = temp;

        if mod(i,ceil(0.001*num_to_run))==0
            temp = toc(t_start);
            fprintf("Finished %d of %d after %s. ETR: %s. ETT: %s\n",i,num_to_run,duration(0,0,temp),duration(0,0,temp/i * (num_to_run-i)),duration(0,0,temp+temp/i * (num_to_run-i)))
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"profiles")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end
else
    %% run in serial
    for i = 1:num_to_run
        current_ind = ind_to_run(i);
        d = D(:,current_ind);
        p = P(:,current_ind);
        profiles(:,current_ind) = profileLikelihood(p,t,d,C,sm,profile_params,profile_likelihood_opts,raw_error_opts);
        if mod(i,ceil(0.001*num_to_run))==0
            temp = toc(t_start);
            fprintf("Finished %d of %d after %s. ETR: %s. ETT: %s\n",i,num_to_run,duration(0,0,temp),duration(0,0,temp/i * (num_to_run-i)),duration(0,0,temp+temp/i * (num_to_run-i)))
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"profiles")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end
end
profiles = reshape(profiles,[npars,cohort_size]);

end
