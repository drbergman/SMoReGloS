function M = finishParameterSetup_Event(M,DT)

if DT < 1e-10
    DT = 0;
end

M.Nsteps = ceil(DT/M.pars.max_dt);
M.dt = DT/M.Nsteps;
