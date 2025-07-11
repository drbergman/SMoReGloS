clearvars;

addpath("../../../src/ProfileLikelihoodFns/")

overwrite_profile = false;
profile_to_clean = "data/Profiles_SMFromABM";

boundary_tolerance = 0.01;
%% load and clean profiles
files.profiles = profile_to_clean;

profiles = cleanProfiles(files);

if overwrite_profile
    file_to_save = profile_to_clean; %#ok<*UNRCH>
else
    file_to_save = profile_to_clean + "_clean.mat";
    copyfile(profile_to_clean + ".mat", file_to_save)
end

save(file_to_save, "profiles", "-append")

%% reset path
rmpath("../../../src/ProfileLikelihoodFns/")
