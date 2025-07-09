function pars = basePars()

pars.max_dt = 4 / 24; % number of days per step

pars.cell_width = 20; % in micrometers; cell is about 20micrometers in diameter
pars.max_tumor_size = Inf; % if tumor exceeds this, stop simulation

%% neighbor parameters
pars.occmax_2d = 6; % below this threshold, a tumor/immune cell can divide; at and above, too many neighbors and so doesn't proliferate
pars.occmax_3d = 20; % below this threshold, a tumor/immune cell can divide; at and above, too many neighbors and so doesn't proliferate

%% tumor parameters
pars.apop_rate = .05; % max death rate of tumor cells (per day)
pars.move_rate_microns = 1 / (1/24); % move rate in microns per day

