function out = sampleFromSMProfiles(x, BS, vals, sm_functional, options)
% runs the SM on a LHS of ODE parameter space as defined by x.
% par_names stores the correspondence between
% indices in x and parameters in the SM. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value. T is
% a dictionary of transformation maps for the parameters. e.g. converting
% to an integer.

n_sm_pars = size(BS,1);

n_abm_pars = length(x);

cm_pars = drawCMParameters(x,options.par_names,options.D,options.T,options.I);

colons = repmat({':'},n_abm_pars,1);
Vq = zeros(n_sm_pars,2);
for pi = 1:n_sm_pars
    Vq(pi,1) = interpn(vals{:},squeeze(BS(pi,colons{:},1)),cm_pars{:});
    Vq(pi,2) = interpn(vals{:},squeeze(BS(pi,colons{:},2)),cm_pars{:});
end

X = lhsdesign(options.nsamps,n_sm_pars,"Smooth","off"); % LHS on [0,1]
points = (1-X).*Vq(:,1)' + X.*Vq(:,2)'; % use these as weights to interpolate between the boundaries

out = zeros(options.nsamps,1);

for i = 1:options.nsamps
    out(i) = sm_functional(points(i,:)');
end

out = options.sum_fn(out);
