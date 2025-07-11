% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 
clear all
clearvars;
addpath("../../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../../src/SurrogateModelFns/")

%% reset persistent variables first
clear rawError solveSM customRawError

%% continue...
make_save = true;

force_serial = false;
n_starts = 1;
assume_independent_time_series = true;

cohort_name = "cohort_2507082020";


[p,lb,ub] = basePars();
npars = numel(p);

% Using base method
% sm.type = "ode";
% sm.solver = @ode45;
% sm.fn = @odefn;
% sm.y0 = [90;10];

% Using custom solve sm method
% sm.custom_solve_sm_fn = @customSolveSM;

% Using custom raw error method
sm.custom_raw_error_fn = @customRawError;

weights = 1;
optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% load ABM data
files.data = sprintf("../../data/%s/summary.mat",cohort_name);

P = optimizeSMParsFromABM(files,sm,p,lb,ub,optim_opts,weights,...
    force_serial=force_serial,n_starts=n_starts,assume_independent_time_series=assume_independent_time_series);

if make_save
    save("data/SMFitToABM.mat","P","cohort_name","sm", "lb", "ub", "optim_opts")
end

%% reset path
rmpath("../../../src/SurrogateModelFns/")
