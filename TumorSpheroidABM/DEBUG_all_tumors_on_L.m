function DEBUG_all_tumors_on_L(M)

assert(all(M.L(M.tumor(:,M.I.ind))~=0))


assert(all(any(M.L(M.tumor(:,M.I.ind))==M.val.tum,2)))