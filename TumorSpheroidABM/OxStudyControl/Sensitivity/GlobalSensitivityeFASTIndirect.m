% A script to run the indirect global sensitivity.

clearvars;

addpath("../../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")
addpath("../../../src/ProfileLikelihoodFns/")
addpath("../../../src/SensitivityFns/")

use_profiles = true;
sort_output = false;

Nr = 15; % number of resamples per factor in ABM space
nsamps = 200; % number of points to sample in LHS for ODE pars
omega_max = 8;
M = 4;
Ns = 65;

files.profiles = "../ProfileLikelihood/data/Profiles_SMFromABM_clean.mat";

D = makeABMParameterDistributionsDictionary();
T = dictionary("occmax_2d",@(x) min(7,floor(x)));

sm_functional = @(p) sum(computeTimeSeries(p, [], [], false, 3));

%% run MOAT
[studied_function,par_names] = setupSampleFromSMFunction(files,sm_functional,T=T,D=D,nsamps=nsamps,use_profiles=use_profiles);
n_abm_pars = length(par_names);
[S1,ST] = efast(studied_function,n_abm_pars,Nr,omega_max,M,Ns);

%% save result
save("data/eFASTIndirect.mat","S1","ST",...
    "par_names","Nr","nsamps","Ns","M","omega_max")

%% clean path
rmpath("../ODEFitting/")
rmpath("../../../src/ProfileLikelihoodFns/")
rmpath("../../../src/SensitivityFns/")


