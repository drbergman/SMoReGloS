% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;
addpath("../../../SurrogateModelFns/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")
addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

%% clear persistent variables
clear rawError solveSM customRawError

file_name = "Profiles_SMFromABM";

files.optimal_parameters = "../ODEFitting/data/SMFitToABM.mat";

load(files.optimal_parameters,"cohort_name")

files.data = sprintf("../../data/%s/summary.mat",cohort_name);
% files.previous_profile_file = "data/Profiles_SMFromABM_clean.mat";

save_all_pars = true;
force_serial = true;
save_every_iter = 50; % wait at least this many iterations between saves
save_every_sec = 5*60; % wait at least this many seconds between saves

load("../ODEFitting/data/SMFitToData_New.mat","sm","lb","ub","optim_opts")
optim_opts.Display = "off";
[p,~,~] = basePars();

%% setup profile params
n_sm_pars = numel(p);

optim_opts.Display = "off";

profile_params.initial_step_prop = .01*ones(n_sm_pars,1); % proportion of best parameter value to begin stepping
profile_params.min_num_steps = 10*ones(n_sm_pars,1); % number of steps to force it to take in the first phase
profile_params.smallest_par_step = [1e-1;1e-1;1e-1]; % do not let the step size go below this as it steps towards the boundary/threshold
profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to the boundary
profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals

profile_params.secondary_step_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after the initial search
profile_params.step_growth_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after successfully extending the profile

% set bounds for optimizing when profiling the other parameters
profile_params.lb = lb;
profile_params.ub = ub;

profile_params.opts = optim_opts;

% specify parameter ranges for bounds on profiling
profile_params.para_ranges = [lb,ub];
profile_params.weights = 1;

%% perform profile
profiles = performProfile(files,sm,profile_params,...
    save_all_pars = save_all_pars, force_serial = force_serial, ...
    save_every_iter = save_every_iter, save_every_sec = save_every_sec );

%% save
if isfield(files,"previous_profile_file")
    files = rmfield(files,"previous_profile_file");
end
save("data/" + file_name,"profiles","files","profile_params")

%% reset path
rmpath("../ODEFitting/")
rmpath("../../../SurrogateModelFns/")
rmpath("../../../ProfileLikelihoodFns/")


