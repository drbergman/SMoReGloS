% A script to run the indirect global sensitivity.

clearvars;

addpath("../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")
addpath("../../src/ProfileLikelihoodFns/")
addpath("../../src/SensitivityFns/")

npoints = 15; % number of points to sample in LHS for ABM pars
nsamps = 100; % number of points to sample in LHS for ODE pars

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

        files.profiles = sprintf("../ProfileLikelihood/data/ProfileLikelihoods_%s_clean.mat",model_type);

        load("../PostAnalysis/data/summary.mat","vals","cohort_size","display_par_names")

        n_abm_pars = length(display_par_names);
        D = makeMOATDistributions(display_par_names);
        T = makeParameterTransformations(display_par_names);

        %% run MOAT
        studied_function = setupSampleFromSMFunction(files,sm_functional,D=D,T=T,nsamps=nsamps,sum_fn=sum_fn,par_names = display_par_names, warnings=false);
        [mu_star,sigma,order] = morris_simple(studied_function,n_abm_pars,npoints);
        display_par_names = display_par_names(order);

        %% save result
        save(sprintf("data/GlobalSensitivityMOATIndirect_%s_%s.mat",model_type,endpoint),"mu_star","sigma","display_par_names","npoints","nsamps")
    end
end
%% clean path
rmpath("../ODEFitting/")
rmpath("../../src/ProfileLikelihoodFns/")
rmpath("../../src/SensitivityFns/")


