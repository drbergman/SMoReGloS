function M = initializeTracked(M)

M.tracked.t = 0;
M.tracked.NT = M.NT;

M.tracked.phases = sum(M.tumor(:,M.I.phase)==1:M.cycle.n_phases,1);

M.tracked.tum_apop = 0;
M.tracked.tum_prolif = 0;
M.tracked.tum_contact_inhibition = 0;
M.tracked.chemo_arrest = zeros(1,M.cycle.n_phases);
M.tracked.arrest_recovery = zeros(1,M.cycle.n_phases);
