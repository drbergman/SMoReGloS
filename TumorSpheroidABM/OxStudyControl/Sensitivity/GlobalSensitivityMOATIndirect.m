clearvars;

% this script will test out how to do sensitivity on the ABM using the SM

%% Program to run

addpath("../../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")
addpath("../ProfileLikelihood/")
addpath("../../../src/SensitivityFns/")
addpath("../../../src/ProfileLikelihoodFns/")

use_profiles = true;
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
%% 2) Parameters : Please fill in
% Maximum number of simulation runs :
% Large number = better estimation of the influence of the factors
% Recommended value : (number of factors + 1) * 10
% The algorithm will maybe exceed this value if it is considered necessary
options.limit_factor = 0.5; % how to set the limit for separating low and high impact factors
options.initialization_factor = 0.8; % how to determine when it has been sufficiently initialized
nsim_max = 1e4;
nsamps = 100; % number of points to sample in LHS for ODE pars
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

% par_names = ["carrying_capacity";"occmax_2d";"move_rate_microns";"g1_to_s";"s_to_g2";"g2_to_m";"m_to_g1"];
files.profiles = "../ProfileLikelihood/data/Profiles_SMFromABM_clean.mat";
D = makeABMParameterDistributionsDictionary();

%% create bounding hypersurfaces
% Number of factors of uncertainty of the function studied :
% num_factors=numel(par_names); 

T = dictionary("occmax_2d", @(x) min(7,floor(x)));
sm_functional = @(p) sum(computeTimeSeries(p,[],[],false,3));
[studied_function,par_names] = setupSampleFromSMFunction(files,sm_functional,T=T,D=D,nsamps=nsamps,use_profiles=use_profiles);
n_abm_pars = length(par_names);
assert(D.numEntries==numel(par_names)) % make sure each parameter has an associated distribution
[mu_star,sigma,order] = morris_simple(studied_function,n_abm_pars,15,sort_output=false);

[~,~,~] = mkdir("data");
save("data/MOATIndirect.mat", "mu_star", "order")

rmpath("../ODEFitting/")
rmpath("../ProfileLikelihood/")
rmpath("../../../src/SensitivityFns/")
rmpath("../../../src/ProfileLikelihoodFns/")