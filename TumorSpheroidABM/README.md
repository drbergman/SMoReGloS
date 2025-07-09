# Agent-based model
The code in this directory implements an extension of a birth-death-migration model of a growing tumor spheroid.
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
If using your own ABM, save the data in an analogous format to plug into the SMoRe GloS framework.
For this example, update the value of `cohort_name = "cohort_2507081432";` to the name of the cohort folder generated in the previous step.

# Fitting the surrogate model
To get optimal surrogate model fits to initialize profile likelihood calculations, run the `OxStudyControl/ODEFitting/FitSMToABM.m` script.
Again, update the `cohort_name = "cohort_2507081432";` to the name of the cohort folder generated in the above step.

# Running the profile likelihood calculations
