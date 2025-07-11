clearvars;

addpath("../../src/myfunctions/")
addpath("../../src/ProfileLikelihoodFns/")
overwrite_profile = false;

model_types = ["exponential", "logistic", "von_bertalanffy"];

for model_type = model_types
    files.profiles = sprintf("data/ProfileLikelihoods_%s",model_type);
    if model_type == "von_bertalanffy"
        boundary_tolerance = [[1e-4;1e-4;1e-4],[1e-1;1;1e-1]];
    else
        boundary_tolerance = 0.0;
    end
    switch model_type
        case "exponential"
            npars = 1;
        case "logistic"
            npars = 2;
        case "von_bertalanffy"
            npars = 3;
    end

    profiles = cleanProfiles(files, boundary_tolerance=boundary_tolerance);

    if overwrite_profile
        file_name = files.profiles; %#ok<*UNRCH>
    else
        file_name = files.profiles + "_clean";
        copyfile(files.profiles + ".mat",file_name + ".mat")
    end

    save(file_name,"profiles","boundary_tolerance","-append")
end
rmpath("../../src/ProfileLikelihoodFns/")

