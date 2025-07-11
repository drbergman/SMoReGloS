% A script to run the indirect global sensitivity.

clearvars;

addpath("../../../src/myfunctions/")
addpath("../../../src/SensitivityFns/")
addpath("../..")

alpha = 0.05; % significance value for CI to determine if enough samples have been computed
ci_relative_spread = 0.1; % how much the confidence interval can spread around the mean of the stochastic simulation output

%% model parameters
Model = allBaseParameters();

Model.setup.ndims = 2;
Model.setup.censor_date = 3;
Model.setup.N0 = 1e2;
Model.setup.agent_initialization_location = "uniform";
Model.setup.carrying_capacity = 1000;

Model.save_pars.make_save = false;
Model.save_pars.dt = Inf;

Model.pars.max_dt = 0.25 / 24; % number of days per step
Model.pars.occmax_3d = 20;
Model.pars.occmax_2d = 5;
Model.pars.apop_rate = 0;
Model.pars.move_rate_microns = 10;

Model.cycle_pars.g1_to_s = 24/11; % * [0.9;1;1.1];
Model.cycle_pars.s_to_g2 = 24/8; % * [0.9;1;1.1];
Model.cycle_pars.g2_to_m = 24/4; % * [0.9;1;1.1];
Model.cycle_pars.m_to_g1 = 24/1; % * [0.9;1;1.1];

Model.chemo_pars.dna_check_g1 = false;
Model.chemo_pars.dna_check_s = false;
Model.chemo_pars.dna_check_g2 = false;
Model.chemo_pars.dna_check_m = false;

Model.chemo_pars.arrest_coeff_g1 = 0.05;
Model.chemo_pars.arrest_coeff_s = 0.05;
Model.chemo_pars.arrest_coeff_g2 = 0.05;
Model.chemo_pars.arrest_coeff_m = 0.05;

Model.plot_pars.plot_fig = false;
Model.plot_pars.plot_location = false;

nsamps = 10;
par_names = ["carrying_capacity";"occmax_2d";"move_rate_microns";"g1_to_s";"s_to_g2";"g2_to_m";"m_to_g1"];
% par_names = ["carrying_capacity";"g1_to_s"];
D = makeABMParameterDistributionsDictionary();

%%

% Nr = 15; % number of resamples per factor in ABM space
Nr = 2; % number of resamples per factor in ABM space
omega_max = 8;
M = 4;
Ns = 65;

n_abm_pars = length(par_names);

%% run eFAST
studied_function = @(x) moatSample(x,Model,par_names,D,nsamps,alpha,ci_relative_spread);
[S1,ST] = efast(studied_function,n_abm_pars,Nr,omega_max,M,Ns);

%% save result
[~,~,~] = mkdir("data");
save("data/eFASTDirect.mat")

%% clean path
rmpath("../../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
rmpath("../../../src/SensitivityFns/")

