function out = sampleFromSMBest(x, OP, vals, sm_functional, options)
% runs the SM on a LHS of ODE parameter space as defined by x.
% par_names stores the correspondence between
% indices in x and parameters in the SM. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value. T is
% a dictionary of transformation maps for the parameters. e.g. converting
% to an integer.

n_sm_pars = size(OP,1);
n_abm_pars = length(x);
cm_pars = drawCMParameters(x,options.par_names,options.D,options.T,options.I);

colons = repmat({':'},n_abm_pars,1);
p_interp = zeros(n_sm_pars,1);
for pi = 1:n_sm_pars
    p_interp(pi) = interpn(vals{:},squeeze(OP(pi,colons{:})),cm_pars{:});
end
out = sm_functional(p_interp);
