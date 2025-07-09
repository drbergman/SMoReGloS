clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

%%
M = allBaseParameters();

%%
M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
% M.setup.carrying_capacity = 1167.726550079491289579891599714756011962890625;
M.setup.carrying_capacity = 500;
M.setup.use_rates_for_intitial_proportions = false;

M.save_pars.make_save = true;
M.save_pars.dt = 1; % set to Inf to not save anything spatial; otherwise dt is in days
M.save_pars.interpolate_tracked = true; % only save certain interpolated values
M.save_pars.t_min = [0 10 24 36 48 72] * 60;
M.save_pars.fields_to_keep = ["t","phases"];

M.pars.max_dt = 0.25 / 24; % number of days per step

% contact inhibition = threshold number of neighboring cells that allow for proliferation
%   neighboring = in Moore neighborhood (8 spots in 2D, 26 in 3D)
%   > threshold ==> cannot proliferate, <= threshold ==> can proliferate (see if phase==M.cycle.m block of performEvents.m)
M.pars.occmax_3d = 20; % contact inhibition when in 3D simulation
M.pars.occmax_2d = 5; % contact inhibition when in 2D simulation
M.pars.move_rate_microns = 8.791732909379968;
M.pars.apop_rate = 0;

M.cycle_pars.g1_to_s = 2.10273160861396490872721187770366668701171875;
M.cycle_pars.s_to_g2 = 2.939904610492844572178228190750814974308013916015625;
M.cycle_pars.g2_to_m = 5.96565977742448705356537175248377025127410888671875;
M.cycle_pars.m_to_g1 = 24.038155802861748355780946440063416957855224609375;

M.flags.arrest_is_death = false;

M.chemo_pars.concentration = 0;

M.chemo_pars.apoptosis_function = "hill";
M.chemo_pars.apop_c0 = 0;
M.chemo_pars.apop_c1 = 1.0;
M.chemo_pars.apop_ec50 = 3;

M.chemo_pars.dna_check_g1 = true;
M.chemo_pars.dna_check_g2 = true;

M.chemo_pars.arrest_function = "hill";

M.chemo_pars.arrest_coeff_g1 = 0.5;
M.chemo_pars.arrest_coeff_g2 = 0.5;
M.chemo_pars.arrest_ec50_g1 = 3;
M.chemo_pars.arrest_ec50_g2 = 3;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = true;
M.plot_pars.make_movie = false;

tic;
M = simPatient(M);
toc;

if isfield(M.save_pars,"sim_identifier")
    fprintf("Finished simulation %s.\n",M.save_pars.sim_identifier)
else
    fprintf("Finished simulation without saving anything.\n")
end
