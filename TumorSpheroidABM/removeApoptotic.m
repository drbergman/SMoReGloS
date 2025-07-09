function M = removeApoptotic(M)

apoptosis_log = M.tumor(:,M.I.event)==2;
apoptosis_ind = find(apoptosis_log);

if ~isempty(apoptosis_ind)
    M.L(M.tumor(apoptosis_ind,M.I.ind)) = 0;
    M.tumor(apoptosis_ind,:)=[]; % get rid of dead cells
end
