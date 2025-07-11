% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;
addpath("../../src/SurrogateModelFns/")
addpath("../../src/ProfileLikelihoodFns/")
addpath("../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")

force_serial = false;
save_all_pars = true;

resample_t = 15:15:75;

model_types = ["exponential", "logistic", "von_bertalanffy"];
files.data = "../PostAnalysis/data/summary.mat";

for model_type = model_types

    %% files
    files.optimal_parameters = sprintf("../ODEFitting/data/OptimalParameters_%s.mat",model_type);
    load(files.optimal_parameters,"sm")

    %% optimization opts
    profile_params.opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    profile_params.weights = 1;

    %% setup SM parameters, bounds, etc
    p = basePars(model_type);
    n_sm_pars = numel(p);

    %% setup profile params
    profile_params.initial_step_prop = .01*ones(n_sm_pars,1);
    profile_params.min_num_steps = 10*ones(n_sm_pars,1);
    profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower boundary
    profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals
    profile_params.secondary_step_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after the initial search
    profile_params.step_growth_factor = 1.5*ones(n_sm_pars,1); % factor by which to increase the step size after successfully extending the profile
    switch model_type
        case "exponential"
            profile_params.smallest_par_step = 0.01; % do not let the step size go below this as it steps towards the boundary/threshold
            profile_params.lb = 0;
            profile_params.ub = 0.2;
            profile_params.para_ranges = [0,1];     % lambda

        case "logistic"
            profile_params.smallest_par_step = [1e-1;50]; % do not let the step size go below this as it steps towards the boundary/threshold

            % set bounds for optimizing when profiling the other parameters
            profile_params.lb = [0;0];
            profile_params.ub = [1;1e6];

            % specify parameter ranges for bounds on profiling
            profile_params.para_ranges = [0,1;     % r
                0,1e6];  % K

        case "von_bertalanffy"
            profile_params.smallest_par_step = [1e-1;1e-1;1e-1]; % do not let the step size go below this as it steps towards the boundary/threshold

            % set bounds for optimizing when profiling the other parameters
            profile_params.lb = [0;1;0];
            profile_params.ub = [100;1e3;100];
            profile_params.A = [-1,0,1]; % -alpha + beta <= 0 <==> alpha >= beta ==> population grows for small enough x
            profile_params.b = 0;

            % specify parameter ranges for bounds on profiling
            profile_params.para_ranges = [0,100;     % alpha
                1,1e3;  % nu
                0,100]; % beta for chemo activating apoptosis
            sm.opts.enforce_inequality = true;
        otherwise
            error("%s is an unspecified SM model.\n",objfn_constants.fn_opts.model_type);
    end

    %% perform profile
    profiles = performProfile(files,sm,profile_params,force_serial=force_serial,...
        save_every_iter=1,save_every_sec=0,...
        save_all_pars=save_all_pars,resample_t=resample_t);

    if isfield(files,"previous_profile_file")
        files = rmfield(files,"previous_profile_file");
    end
    save(sprintf("data/ProfileLikelihoods_%s.mat",model_type),"profiles","files","profile_params")
end
%% reset path
rmpath("../ODEFitting/")
rmpath("../../src/SurrogateModelFns/")
rmpath("../../src/ProfileLikelihoodFns/")
