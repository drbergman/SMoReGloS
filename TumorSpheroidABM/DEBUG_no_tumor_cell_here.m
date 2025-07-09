function DEBUG_no_tumor_cell_here(M,ind)

assert(all(M.tumor(:,M.I.ind)~=ind))