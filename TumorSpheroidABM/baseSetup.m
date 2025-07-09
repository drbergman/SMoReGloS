function setup = baseSetup()

% set base values for setup parameters

setup.use_carrying_capacity_for_grid_size = true;
setup.carrying_capacity = 6500;

setup.grid_size_microns_x = 400;
setup.grid_size_microns_y = 400;
setup.grid_size_microns_z = 400;

setup.ndims = 3;

setup.censor_date = 25;
setup.N0 = 1000;

setup.agent_initialization_location = "uniform";

setup.c = -3.603357085551339; % see Test_findInitialization for code that got to this number
setup.e = 2.986939791722032; % see Test_findInitialization for code that got to this number

setup.use_rates_for_intitial_proportions = true; % whether to use transition rates to approximate initial phase distributions or to use a fixed value based on data
