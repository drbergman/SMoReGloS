% script that calls the MOAT routines for sensitivity of the ABM parameters

%% Program to run

addpath("../../../src/myfunctions/")
addpath("../..")

%
% This algorithm is an adaptation of the method of Sensitivity Analysis
% called the Morris method.
%
% Sensitivity analysis is used to estimate the influence of uncertainty
% factors on the output of a function.
% The Morris method is sometimes referenced to as a qualitative method : it
% gives rough estimations with a limited number of calculations.
% The Morris method can be used to simplify a function, as a first step. It
% can identify the factors with a low influence which can be fixed.
% For further information :
% Saltelli, A., Tarantola, S., Campolongo, F., and Ratto, M. (2004). Sensitivity Analysis in Practice - A Guide to Assessing Scientific Models. Wiley.
% 
% This algorithm reduces the risk to underestimate and fix non-negligible
% factors. It is presented in:
% Henri Sohier, Helene Piet-Lahanier, Jean-Loup Farges, Analysis and optimization of an air-launch-to-orbit separation, Acta Astronautica, Volume 108, March泡pril 2015, Pages 18-29, ISSN 0094-5765, http://dx.doi.org/10.1016/j.actaastro.2014.11.043.
% (http://www.sciencedirect.com/science/article/pii/S0094576514004974)
%
% This program is divided in 6 parts:
% 1) Clearing the memory
% 2) Parameters : Please fill in
% 3) Initialization of the variables
% 4) Loop
% 5) Output text
% 6) Figure
%
% Please fill in the second part to apply the algorithm to your function.
% Do not change the parameters to see the results with the "modified Sobol
% test function".
%
% This program outputs a figure as well as a short summary in the console.
% Consider fixing the factors which appear as negligible on the left of the
% figure (not necessary all the factors under the limit).
%% 1) Clearing the memory
clearvars; % Clears the memory

%% 2) Parameters : Please fill in
% Maximum number of simulation runs :
% Large number = better estimation of the influence of the factors
% Recommended value : (number of factors + 1) * 10
% The algorithm will maybe exceed this value if it is considered necessary
alpha = 0.05; % significance value for CI to determine if enough samples have been computed
options.limit_factor = 0.5; % how to set the limit for separating low and high impact factors
options.initialization_factor = 0.8; % how to determine when it has been sufficiently initialized
ci_relative_spread = 0.1; % how much the confidence interval can spread around the mean of the stochastic simulation output
nsim_max = 210;
% Function studied :
% Replace test_function by the name of your function. It must be a 
% function with one multidimensional input x. x must represent the values 
% of the uncertainty factors in the quantiles hyperspace (the i-th 
% coordinate of x is not the actual value of the i-th factor, but the 
% corresponding value of the cumulative distribution function of the i-th 
% factor). To adapt your function, first calculate the actual values of 
% the factors by applying the inverse of their cumulative distribution 
% function to each coordinate of x; Matlab includes such inverses: 
% mathworks.com/help/stats/icdf.html ).
M = allBaseParameters();

M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
M.setup.carrying_capacity = 1000;

M.save_pars.make_save = false;
M.save_pars.dt = Inf;

M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.occmax_3d = 20;
M.pars.occmax_2d = 5;
M.pars.apop_rate = 0;
M.pars.move_rate_microns = 10;

M.cycle_pars.g1_to_s = 24/11; % * [0.9;1;1.1];
M.cycle_pars.s_to_g2 = 24/8; % * [0.9;1;1.1];
M.cycle_pars.g2_to_m = 24/4; % * [0.9;1;1.1];
M.cycle_pars.m_to_g1 = 24/1; % * [0.9;1;1.1];

M.chemo_pars.dna_check_g1 = false;
M.chemo_pars.dna_check_s = false;
M.chemo_pars.dna_check_g2 = false;
M.chemo_pars.dna_check_m = false;

M.chemo_pars.arrest_coeff_g1 = 0.05;
M.chemo_pars.arrest_coeff_s = 0.05;
M.chemo_pars.arrest_coeff_g2 = 0.05;
M.chemo_pars.arrest_coeff_m = 0.05;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

nsamps = 10;
par_names = ["carrying_capacity";"occmax_2d";"move_rate_microns";"g1_to_s";"s_to_g2";"g2_to_m";"m_to_g1"];
% par_names = ["carrying_capacity";"g1_to_s"];
D = makeABMParameterDistributionsDictionary();

% Number of factors of uncertainty of the function studied :
nfac=numel(par_names); 

assert(nfac==numel(par_names)) % make sure that there is a value for each of the parameters to be varied
assert(D.numEntries==numel(par_names)) % make sure each parameter has an associated distribution
studied_function = @(x) moatSample(x,M,par_names,D,nsamps,alpha,ci_relative_spread);
[mu_star,sigma,order] = morris_simple(studied_function,nfac,15);

mkdir("data")
save("data/MOATDirect.mat", "mu_star", "order", "par_names")

rmpath("../..")
