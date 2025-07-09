function DEBUG_all_tum_L_in_tumors(M)

L_ind = find(any(M.L==M.val.tum,2));
T_ind = sort(M.tumor(:,M.I.ind));
assert(isempty(setdiff(L_ind,T_ind)))