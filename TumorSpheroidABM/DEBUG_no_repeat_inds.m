function DEBUG_no_repeat_inds(M)

assert(length(unique(M.tumor(:,M.I.ind)))==size(M.tumor,1))