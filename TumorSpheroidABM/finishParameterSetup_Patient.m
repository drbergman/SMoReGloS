function M = finishParameterSetup_Patient(M)

if M.setup.use_carrying_capacity_for_grid_size
    K = M.setup.carrying_capacity * M.pars.cell_width^M.setup.ndims; % carrying capacity in cubic microns
    M.setup.grid_size_microns_x = ceil(nthroot(K,M.setup.ndims));
    M.setup.grid_size_microns_y = ceil(nthroot(K/M.setup.grid_size_microns_x,M.setup.ndims-1));
    if M.setup.ndims==3
        M.setup.grid_size_microns_z = ceil(K/(M.setup.grid_size_microns_x*M.setup.grid_size_microns_y));
    end
end

if M.setup.ndims == 3
    M.setup.grid_size_microns = [M.setup.grid_size_microns_x,M.setup.grid_size_microns_y,M.setup.grid_size_microns_z];

    %% neighbor stuff
    M.pars.neighbors = allCombos(-1:1,-1:1,-1:1,'matlab');
    M.pars.neighbors(all(M.pars.neighbors==0,2),:) = []; % don't count self as neighbor
    M.pars.neighbors_VN = M.pars.neighbors;
    M.pars.neighbors_VN(sum(abs(M.pars.neighbors_VN),2)>1,:) = [];

    M.pars.neighbor_weights = 1./sqrt(sum(M.pars.neighbors.^2,2));
    M.pars.left_neighbors = M.pars.neighbors(:,1) == -1;
    M.pars.right_neighbors = M.pars.neighbors(:,1) == 1;
    M.pars.front_neighbors = M.pars.neighbors(:,2) == -1;
    M.pars.back_neighbors = M.pars.neighbors(:,2) == 1;
    M.pars.bottom_neighbors = M.pars.neighbors(:,3) == -1;
    M.pars.top_neighbors = M.pars.neighbors(:,3) == 1;

    M.pars.left_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == -1);
    M.pars.right_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == 1);
    M.pars.front_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == -1);
    M.pars.back_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == 1);
    M.pars.bottom_neighbors_VN_ind = find(M.pars.neighbors_VN(:,3) == -1);
    M.pars.top_neighbors_VN_ind = find(M.pars.neighbors_VN(:,3) == 1);

    M.pars.occmax = M.pars.occmax_3d;

else
    M.setup.grid_size_microns = [M.setup.grid_size_microns_x,M.setup.grid_size_microns_y];

    %% neighbor stuff
    M.pars.neighbors = allCombos(-1:1,-1:1,'matlab');
    M.pars.neighbors(all(M.pars.neighbors==0,2),:) = []; % don't count self as neighbor
    M.pars.neighbors_VN = M.pars.neighbors;
    M.pars.neighbors_VN(sum(abs(M.pars.neighbors_VN),2)>1,:) = [];

    M.pars.neighbor_weights = 1./sqrt(sum(M.pars.neighbors.^2,2));
    M.pars.left_neighbors = M.pars.neighbors(:,1) == -1;
    M.pars.right_neighbors = M.pars.neighbors(:,1) == 1;
    M.pars.front_neighbors = M.pars.neighbors(:,2) == -1;
    M.pars.back_neighbors = M.pars.neighbors(:,2) == 1;

    M.pars.left_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == -1);
    M.pars.right_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == 1);
    M.pars.front_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == -1);
    M.pars.back_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == 1);

    M.pars.occmax = M.pars.occmax_2d;

end

M.pars.move_rate = M.pars.move_rate_microns / M.pars.cell_width;

%% index and value stuff
M.I = buildIndices(M.setup.ndims);
M.val = buildLatticeVals();
M.cycle = buildCycle();

%% cycle stuff
M.cycle.phase_names = string(fieldnames(M.cycle));
M.cycle.n_phases = length(M.cycle.phase_names); % number of phases included in current model

M.chemo_pars.dna_check = false(1,M.cycle.n_phases);
M.chemo_pars.arrest_prob = zeros(1,M.cycle.n_phases);
for i = 1:M.cycle.n_phases
    M.chemo_pars.dna_check(M.cycle.(M.cycle.phase_names(i))) = M.chemo_pars.(sprintf("dna_check_%s",M.cycle.phase_names(i)));
    switch M.chemo_pars.arrest_function
        case "linear"
            M.chemo_pars.arrest_prob(M.cycle.(M.cycle.phase_names(i))) = M.chemo_pars.(sprintf("arrest_coeff_%s",M.cycle.phase_names(i))) * M.chemo_pars.concentration;
        case "hill"
            M.chemo_pars.arrest_prob(M.cycle.(M.cycle.phase_names(i))) = M.chemo_pars.(sprintf("arrest_coeff_%s",M.cycle.phase_names(i))) * M.chemo_pars.concentration / (M.chemo_pars.("arrest_ec50_" + M.cycle.phase_names(i)) + M.chemo_pars.concentration);
        otherwise
            error("%s is not a specified arrest function.\n",M.chemo_pars.arrest_function)
    end
end

M.chemo_pars.dna_check = M.chemo_pars.dna_check & M.chemo_pars.arrest_prob > 0; % in case the check is on, but the probability is 0; will save on drawing random numbers in these cases

M.cycle.advancer = [M.cycle.s,M.cycle.g2,M.cycle.m,M.cycle.g1,M.cycle.g1,M.cycle.g1,M.cycle.g1,M.cycle.g1]; % this gives the index of the next cycle (final four are for arrested compartments returning to g1)
M.cycle_pars.transition_rates = [M.cycle_pars.g1_to_s,M.cycle_pars.s_to_g2,M.cycle_pars.g2_to_m,M.cycle_pars.m_to_g1,M.cycle_pars.g1a_to_g1,M.cycle_pars.sa_to_g1,M.cycle_pars.g2a_to_g1,M.cycle_pars.ma_to_g1];

%% apoptosis rates
switch M.chemo_pars.apoptosis_function
    case "constant" % arrested cells have apoptosis rate that is independent of the chemo concentration
        M.chemo_pars.apop_rate = M.chemo_pars.apop_c0;
    case "linear" % arrested cells have apoptosis rate that grows linearly with the chemo concentration
        M.chemo_pars.apop_rate = M.chemo_pars.apop_c0 + M.chemo_pars.apop_c1*M.chemo_pars.concentration;
    case "hill" % arrested cells have apoptosis rate that grows like a hill function with the chemo concentration
        M.chemo_pars.apop_rate = M.chemo_pars.apop_c0 + M.chemo_pars.apop_c1*M.chemo_pars.concentration/(M.chemo_pars.apop_ec50+M.chemo_pars.concentration);
    otherwise
        error("%s is not a specified apoptosis function.\n",M.chemo_pars.apoptosis_function)
end
M.flags.cycling_apoptosis_is_nonzero = M.pars.apop_rate>0; % whether apoptosis is possible for cycling cells
M.flags.arrested_apoptosis_is_nonzero = M.chemo_pars.apop_rate>0; % whether apoptosis is possible for arrested cells

%% events
M = initializeEvents(M);
