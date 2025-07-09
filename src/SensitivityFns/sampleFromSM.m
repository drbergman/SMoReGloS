% function out = sampleFromSM(x,par_names,BS,T,D,vals,nsamps,sm_functional,C,fn_opts,sum_fn)
function out = sampleFromSM(x, files, vals, sm_functional, options)
% runs the SM on a LHS of ODE parameter space as defined by x.
% par_names stores the correspondence between
% indices in x and parameters in the SM. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value. T is
% a dictionary of transformation maps for the parameters. e.g. converting
% to an integer.

arguments
    x (:,1) double
    files struct % structure containing the files with the necessary info for this
    % BS double % at some point, change this so the user passes in the file that contains this info and build the BS from there
    vals cell % at some point, change this so the user passes in the file that contains this info and build the vals from there
    sm_functional function_handle
    options.par_names {mustBeText} = strings(1,numel(x))
    options.T dictionary = configureDictionary("string","function_handle")
    options.D dictionary = configureDictionary("string","prob.ProbabilityDistribution")
    options.nsamps double {mustBeInteger} = 100
    options.sum_fn function_handle = @mean
    options.use_profiles logical = true
end
persistent BS;
persistent OP;
if options.use_profiles
    if isempty(BS)
        BS = loadBoundingSurfaces(files.profiles);
    end
    n_sm_pars = size(BS,1);
else
    if isempty(OP)
        OP = loadOptimalParameters(files);
    end
    n_sm_pars = size(OP,1);
end

n_abm_pars = length(x);
abm_pars = cell(n_abm_pars,1);
for i = 1:n_abm_pars
    if ~any(options.par_names(i) == options.D.keys)
        continue;
    end
    abm_pars{i} = icdf(options.D(options.par_names(i)),x(i));
    if any(options.par_names(i) == options.T.keys)
        abm_pars{i} = feval(options.T(options.par_names(i)),abm_pars{i});
    end
end

colons = repmat({':'},n_abm_pars,1);
if options.use_profiles
    Vq = zeros(n_sm_pars,2);
    for pi = 1:n_sm_pars
        Vq(pi,1) = interpn(vals{:},squeeze(BS(pi,colons{:},1)),abm_pars{:});
        Vq(pi,2) = interpn(vals{:},squeeze(BS(pi,colons{:},2)),abm_pars{:});
    end

    X = lhsdesign(options.nsamps,n_sm_pars,"Smooth","off"); % LHS on [0,1]
    points = (1-X).*Vq(:,1)' + X.*Vq(:,2)'; % use these as weights to interpolate between the boundaries

    out = zeros(options.nsamps,1);

    for i = 1:options.nsamps
        out(i) = sm_functional(points(i,:)');
    end

    out = options.sum_fn(out);
else
    p_interp = zeros(n_sm_pars,1);
    for pi = 1:n_sm_pars
        p_interp(pi) = interpn(vals{:},squeeze(OP(pi,colons{:})),abm_pars{:});
    end
    out = sm_functional(p_interp);

end

end

function P = loadOptimalParameters(files)

if ~isfield(files,"optimal_parameters")
    load(files.profiles,"files")
end
load(files.optimal_parameters,"P")

end

function BS = loadBoundingSurfaces(profile_file)

load(profile_file,"profiles","files");
n_sm_pars = size(profiles,1);
profiles = reshape(profiles,n_sm_pars,[]);
n_abm_vecs = size(profiles,2);


BS = zeros(n_sm_pars,n_abm_vecs,2);
threshold = chi2inv(0.95,n_sm_pars);
for i = 1:n_abm_vecs
    for j = 1:n_sm_pars
        [BS(j,i,1),BS(j,i,2)] = getProfileBounds(profiles{j,i}([j,end],:),threshold);
    end
end

load(files.data,"cohort_size")
BS = reshape(BS,[n_sm_pars,cohort_size,2]);

end