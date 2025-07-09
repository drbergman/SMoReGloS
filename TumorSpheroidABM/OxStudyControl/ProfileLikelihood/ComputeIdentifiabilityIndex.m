clearvars;

addpath("../../../ProfileLikelihoodFns/")


files.profiles = sprintf("data/Profiles_SMFromABM_clean.mat");
indices = identifiabilityIndex(files);

file_name = sprintf("data/IdentifiabilityIndex");
if exist(file_name,"file") && ~overwrite_save_files
    return
end
save(file_name,"files","indices")

rmpath("../../../ProfileLikelihoodFns/")
