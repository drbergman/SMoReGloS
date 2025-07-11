# Agent-based model
The code in this directory implements the 3D vascular growth model.
The subdirectories carry out the steps required to run the SMoRe GloS framework on this model.
**For each step, make sure the directory containing the script is the current working directory in MATLAB.**

# Running the initial simulations
The `src/RunCohort.m` function runs all the simulations necessary to construct the hypersurfaces for SMoRe GloS.
These 486 simulations (3x3x3x3=81 ABM parameter vectors with 6 replicates each) take a considerable amount of time to run (the authors have not run them all locally).
The `src/runSamples.sbat` sbatch script can be used on a SLURM-based HPC with appropriate modifications for username, account, etc.

# Summarizing the results
The scripts `PostAnalysis/CollectSimData.m` and `PostAnalysis/SummarizeSimData.m` should be run consecutively to process the data into the necessary format for downstream analysis in the SMoRe GloS framework.

# Fitting the surrogate model
To get optimal surrogate model fits to initialize profile likelihood calculations, run the `ODEFitting/FitODEToABM.m` script.
This will compute them for all three surrogate models.

# Running the profile likelihood calculations
To compute the profile likelihoods, run the `ProfileLikelihood/ProfileSMFromABM.m` script.
This will compute them for all three surrogate models and will take about 5 minutes.

# Clean the profiles
To restrict the profiles to the 95% confidence intervals, run `ProfileLikelihood/CleanProfiles.m`.

# Compute identifiability index
The script `ProfileLikelihood/ComputeIdentifiabilityIndex.m` computes the identifiability indices for each SM parameter at all the sampled ABM parameter vectors.

# Compute sensitivity
To compute GSA directly, the steps to create the ABM parameter sampling and running the simulations is split in two pieces: first use `Sensitivity/MakeMOATSamplePars.m` and `Sensitivity/MakeEFASTSamplePars.m` to sample the ABM parameter vectors, and then use `Sensitivity/RunMOATCohort.m` and `Sensitivity/RunEFASTCohort.m` to run the ABM simulations locally.

Completing these cohorts will take significant compute time when run locally.
Sbatch scripts are provided that can help set up and run the simulations on an HPC: `Sensitivity/runMOATSamples.sbat` and `Sensitivity/runEFASTSamples.sbat`.

The `Sensitivity/GlobalSensitivityMOATDirect.m` and `Sensitivity/GlobalSensitivityEFASTDirect.m` scripts collate these simulations and perform the GSA.

To use SMoRe GloS, run the `Sensitivity/GlobalSensitivityMOATIndirect.m` and `Sensitivity/GlobalSensitivityEFASTIndirect.m` scripts.
These compute the GSA using all three surrogate models.

Finally, the `Sensitivity/CompareAllMOAT.m` and `Sensitivity/CompareAllEFAST.m` scripts create several plots comparing the two methods for GSA across the three endpoints.