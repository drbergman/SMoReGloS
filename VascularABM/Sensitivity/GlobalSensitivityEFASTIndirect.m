% A script to run the indirect global sensitivity.

clearvars;

addpath("../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")
addpath("../../src/ProfileLikelihoodFns/")
addpath("../../src/SensitivityFns/")

Nr = 15; % number of resamples per factor in ABM space
nsamps = 200; % number of points to sample in LHS for ODE pars
omega_max = 8;
M = 4;
Ns = 249;

sort_output = false;

model_types = ["exponential", "logistic", "von_bertalanffy"];
endpoints = ["final_size", "AUC", "time_to_half"];

[~,~,~] = mkdir("data");

for model_type = model_types
    for endpoint = endpoints

        if model_type == "von_bertalanffy" && endpoint == "time_to_half"
            continue;
        end

        switch endpoint
            case "final_size"
                sm_functional = @(p) computeSMEndpoint(p,[],model_type);
            case "AUC"
                sm_functional = @(p) computeSMAUC(p,[],model_type);
            case "time_to_half"
                sm_functional = @(p) computeTimeToHalf(p,model_type);
        end

        switch model_type
            case "exponential"
                sum_fn = @mean;
            case "logistic"
                sum_fn = @mean; % use this for logistic growth
            case "von_bertalanffy"
                if endpoint == "final_size"
                    sum_fn = @summarizeSMLHS; % use this for the VB model when the parameters make it likely that simulations blow up
                elseif endpoint == "AUC"
                    sum_fn = @mean;
                end
        end

        files.profiles = sprintf("../ProfileLikelihood/data/ProfileLikelihoods_%s.mat",model_type);
        files.data = "../PostAnalysis/data/summary.mat";

        load("../PostAnalysis/data/summary.mat","display_par_names")

        display_par_names_original = display_par_names;
        n_abm_pars = length(display_par_names_original);
        D = makeMOATDistributions(display_par_names_original);
        T = makeParameterTransformations(display_par_names_original);

        %% run eFAST
        studied_function = setupSampleFromSMFunction(files, sm_functional, D=D,T=T,nsamps=nsamps,par_names=display_par_names,sum_fn=sum_fn,warnings=false);
        [S1,ST] = efast(studied_function,n_abm_pars,Nr,omega_max,M,Ns);

        %% save result
        save(sprintf("data/GlobalSensitivityEFASTIndirect_%s_%s.mat",model_type,endpoint),"S1","ST","display_par_names","Nr","nsamps","Ns","M","omega_max")

    end
end
%% clean path
rmpath("../ODEFitting/")
rmpath("../../src/ProfileLikelihoodFns/")
rmpath("../../src/SensitivityFns/")


