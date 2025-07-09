function M = saveInitialModelData(M)

sim_identifier = M.save_pars.sim_identifier;
save(sprintf("data/sims/%s/output_constants",M.save_pars.sim_identifier),'sim_identifier','-v7.3') % start saving this; as new cohorts add this sim, the sim_ids from that will be appended here
clear sim_identifier

if M.save_pars.save_constants
    grid_size = M.grid.size;
    setup = M.setup;
    pars = M.pars;
    cycle_pars = M.cycle_pars;
    chemo_pars = M.chemo_pars;
    plot_pars = M.plot_pars;
    flags = M.flags;
    save_pars = M.save_pars;

    save(sprintf("data/sims/%s/output_constants",M.save_pars.sim_identifier),'-regexp','[^M]','-append')
end


if M.save_pars.dt < Inf
    M = saveModelData(M);
end
