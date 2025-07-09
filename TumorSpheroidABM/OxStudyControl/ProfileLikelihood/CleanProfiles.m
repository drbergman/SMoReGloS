clearvars;

addpath("../../../ProfileLikelihoodFns/")

overwrite_profile = false;
files.profiles = "data/Profiles_SMFromABM_clean";

boundary_tolerance = 0.01;
%% load and clean profiles

profiles = cleanProfiles(files,boundary_tolerance=boundary_tolerance);

if overwrite_profile
    save(profile_to_clean,"profiles") %#ok<*UNRCH>
else
    save(profile_to_clean + "_clean","profiles")
end

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
