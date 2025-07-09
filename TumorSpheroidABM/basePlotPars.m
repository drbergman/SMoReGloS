function plot_pars = basePlotPars()

%% Plotting properties

plot_pars.plot_fig = false;
plot_pars.plot_location = false; % don't bother plotting where all the cells are. it takes a long time
plot_pars.make_movie = false;
plot_pars.plot_every = 1;
plot_pars.plot_offset = 0;

% each axis: [-m m]
plot_pars.m=30;

plot_pars.default_font_size = 16;
