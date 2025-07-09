function LP = linkByName(LP,cohort_pars)

for li = 1:numel(cohort_pars.linkings) % linking index
    index_in_LP = [];
    vals = {};
    paths = {};
    for i = 1:numel(LP)
        if ~iscell(LP(i).path) && any(LP(i).path(end) == cohort_pars.linkings{li})
            index_in_LP(end+1) = i;
            vals{end+1} = LP(i).values;
            paths{end+1} = LP(i).path;
        end
    end
    if length(index_in_LP)>1 % only make this change if two arrest coeffs are varied
        assert(numel(unique(cellfun(@numel,vals)))==1) % make sure that each of these values has the same number of elements
        LP(end+1).values = cat(2,vals{:});
        LP(end).path = paths;
        LP(index_in_LP) = [];
    end
end
