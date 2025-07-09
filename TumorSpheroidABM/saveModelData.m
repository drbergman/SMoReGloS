function M = saveModelData(M)

time = M.t;

tumor_locations = M.save_pars.integrify(M.tumor(:,M.I.subs));
tumor_phases = uint8(M.tumor(:,M.I.phase));

save(sprintf("data/sims/%s/output_%08d",M.save_pars.sim_identifier,M.save_index),'-regexp','[^M]','-v7.3')

M.next_save_time = M.next_save_time + M.save_pars.dt;
M.save_index = M.save_index + 1;