function T = makeParameterTransformations(par_names)

% creates a dictionary of parameter transformations. most will be the
% identity. some will convert to integers. Others may do other things yet
% unknown.

T = dictionary();

for i = 1:numel(par_names)
    switch par_names(i)
        case "p_{lim}"
            T("p_{lim}") = @(x) round(x);
        otherwise
            T(par_names(i)) = @(x) x;
    end
end