clearvars;

% This script will go through the sims in data/sims and make sure they all
% belong to at least one cohort. If not, they will be deleted

delete_sims = true; % whether to actually delete these sims, or just identify them

f = dir("data/sims/*");
f = f([f.isdir]);

for i = numel(f):-1:1
    if ~startsWith(f(i).name,digitsPattern(1)) % I only want folders that start with a number (not "." or ".." etc)
        f(i) = [];
    else
        f(i).name = string(f(i).name);
    end
end

ids = [f.name];

cohorts = dir("data/cohort_*");

for i = 1:numel(cohorts)
    C = load(sprintf("data/%s/output.mat",cohorts(i).name),"ids");
    ids = setdiff(ids,C.ids);
    if isempty(ids)
        break;
    end
end

if delete_sims
    for i = 1:numel(ids)
        rmdir(sprintf("data/sims/%s",ids(i)),'s');
    end
else
    if isempty(ids)
        fprintf("All simulations belong to a cohort. That's good.\n");
    else
        fprintf("Found %d simulations that do not belong to a cohort. These may be disposable.\n",numel(ids))
    end
end
