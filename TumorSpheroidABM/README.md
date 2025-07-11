# Agent-based model
The code in this directory implements the 2D in vitro proliferation assay ABM.
The subdirectories carry out the steps required to run the SMoRe GloS framework on this model.
**For each step, make sure the directory containing the script is the current working directory in MATLAB.**

# Running the initial simulations
The `RunCohort.m` script runs the initial simulations of the model used to set up the mapping between the ABM and the surrogate model.
For each of 7 parameters, 3 values are chosen of each are used to run the model for a total of 2187 different parameter combinations.
At each parameter combination, 6 simulations are run, each with a different random seed.
This produces a total of 13122 simulations.
The simulations take around 7 minutes to complete on a 2024 MacBook Pro with an Apple M4 Max chip.

Note: `cohort_pars.min_parfor_num = Inf;` can be set to a lower value to run the simulations in parallel if the Parallel Computing Toolbox is available.

At the end of the simulations, the console will display the name of the cohort folder generated, which will be in `TumorSpheroidABM/data/` and will be in the format `cohort_YYMMDDHHMM` where `YYMMDDHHMM` is the date and time at which the cohort was started.

# Summarizing the results
The `OxStudyControl/PostAnalysis/SummarizeCohort.m` script will summarize the results of the ABM and write them to a file that can be used by the rest of the SMoRe GloS framework.
This file is located within the cohort folder at `TumorSpheroidABM/data/cohort_YYMMDDHHMM/summary.mat`.
If using your own ABM, save the data in an analogous format to plug into the SMoRe GloS framework.
For this example, update the value of `cohort_name = "cohort_2507081432";` to the name of the cohort folder generated in the previous step.

# Fitting the surrogate model
To get optimal surrogate model fits to initialize profile likelihood calculations, run the `OxStudyControl/ODEFitting/FitSMToABM.m` script.
Again, update the `cohort_name = "cohort_2507081432";` to the name of the cohort folder generated in the above step.
The optimal fits will be saved at `TumorSpheroidABM/OxStudyControl/ODEFitting/data/SMFitToABM.mat`.

# Running the profile likelihood calculations
To compute the profile likelihoods, run the `OxStudyControl/ProfileLikelihood/ProfileSMFromABM.m` script.
This script assumes that the output of the previous step is left in place, so no changes need to be made, i.e., there is no `cohort_name` variable to set.
This step takes approximately 20 minutes on a MacBook Pro (13-inch, M1, 2020) equipped with an Apple M1 chip featuring an 8-core CPU and 8 GB of unified memory.
Intermediate saves are made regularly throughout.
If it fails partway through, set the `% files.previous_profile_file = "data/temp_profile_01.mat";` variable to match the temporary file saved in `ProfileLikelihood/data`.
This script writes out the profiles to `ProfileLikelihood/data/Profiles_SMFromABM.mat`.

# Clean the profiles
The profiles as created by the `performProfile.m` function called in `ProfileSMFromABM.m` explore 1D parameter space until the parameter either exceeds the user-specified bounds or the chi-squared value exceeds the threshold value.
In the latter case, this function does not determine the parameter value at which the threshold is reached, instead just recording the parameter and chi-squared values.
The downstream steps of SMoRe GloS relies on these intervals being accurate and this step "cleans" the profiles so that they actually represent the 95% confidence intervals.

To do so, run `OxStudyControl/ProfileLikelihood/CleanProfiles.m`.
Make sure the `profile_to_clean` is set to the name of the file created in the previous step (only necessary if the file was renamed).
This script is set up to save the cleaned profiles to a new file with `_clean` appended to the filename (before the extension).

# Compute identifiability index
The script `OxStudyControl/ProfileLikelihood/ComputeIdentifiabilityIndex.m` computes the identifiability indices for each SM parameter at all the sampled ABM parameter vectors.
Compute the percentage of each value using calls such as

```
for i = 1:3 % index values
    mean(indices == i, 2); % ABM parameter vectors correspond to columns
end
```

# Compute sensitivity
The scripts `OxStudyControl/Sensitivity/GlobalSensitivityMOATDirect.m`, etc. can be run to perform GSA.
The `CompareMOAT.m` and `CompareEFAST.m` scripts in the same folder compare the direct and indirect methods.