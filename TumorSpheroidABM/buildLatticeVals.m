function val = buildLatticeVals()

%% Set up occupancy array
val.tum = 1; % vals to be put in L for tumor cells
val.tum_apop = 2; % val to be put in L for apoptotic tumor cells
val.unused = 5; % val to be put in L for free spaces in computeQ for indexing purposes
