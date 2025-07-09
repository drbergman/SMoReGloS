function [out,par_names] = setupSampleFromSMFunction(files, sm_functional, options)

arguments
    files struct % structure containing the files with the necessary info for this
    sm_functional function_handle
    options.warnings = true;
    options.par_names {mustBeText} = strings(0,1)
    options.D dictionary = configureDictionary("string","prob.ProbabilityDistribution")
    options.T dictionary = configureDictionary("string","function_handle")
    options.I dictionary = configureDictionary("string","double") % index of parameters
    options.nsamps double {mustBeInteger} = 100
    options.sum_fn function_handle = @mean
    options.use_profiles logical = true
end

vals = loadABMValues(files);
nfacs = numel(vals);
if isempty(options.par_names) && options.D.numEntries~=0
    % then expect that parameter names are needed, not explicitly provided, and go find them
    if isfield(files,"data")
        temp = load(files.data,"par_names");
    else
        temp = load(files.profiles,"files");
        temp = load(temp.files.data,"par_names");
    end
    options.par_names = temp.par_names;
elseif options.warnings
    warning_message = sprintf("User supplied parameter names.\n" + ...
            "\tMake sure they are in the same order as they are varied in the CM cohort.\n" + ...
            "\tDisable this warning by setting setupSampleFromSMFunction(..., warnings=false)");
    warning(warning_message) % warning doesn't do `\n` and `\t` correctly though...
end

par_names = options.par_names;
if ~isempty(options.par_names) && options.I.numEntries<numel(options.par_names)
    % set up index dictionary
    for i = 1:length(options.par_names)
        options.I(options.par_names(i)) = i; % the parameters in par_names are expected to be in the order that they vary in the cohort dimensions
    end
end
    

if options.use_profiles
    BS = loadBoundingSurfaces(files.profiles);
    out = @(x) sampleFromSMProfiles(x, BS, vals, sm_functional, options);

else
    OP = loadOptimalParameters(files);
    out = @(x) sampleFromSMBest(x, OP, vals, sm_functional, options);
end
end

function vals = loadABMValues(files)

if ~isfield(files,"data")
    load(files.profiles,"files")
end
load(files.data,"vals")

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

function P = loadOptimalParameters(files)

if ~isfield(files,"optimal_parameters")
    load(files.profiles,"files")
end
load(files.optimal_parameters,"P")

end