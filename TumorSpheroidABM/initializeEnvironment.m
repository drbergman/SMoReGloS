function M = initializeEnvironment(M)

M = initializeGrid(M);

M = initializeTumor(M);

M.L = zeros(M.V_tot,1);
M.L(M.tumor(:,M.I.ind)) = M.val.tum;
