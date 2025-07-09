function cm_pars = drawCMParameters(x,par_names,D,T,I)

cm_pars = cell(length(x),1);
for i = 1:length(x)
    if ~any(par_names(i) == D.keys)
        continue;
    end
    idx = I(par_names(i));
    cm_pars{idx} = icdf(D(par_names(i)),x(i));
    if any(par_names(i) == T.keys)
        cm_pars{idx} = feval(T(par_names(i)),cm_pars{idx});
    end
end

end