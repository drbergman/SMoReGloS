function D = makeMOATDistributions(par_names)

% defines distributions on each of the varied abm parameters

D = dictionary();
for i = 1:numel(par_names)
    switch par_names(i)
        case "p_{div}"
            D("p_{div}") = makedist("Uniform",0.05,0.245);
            dmin = 0; dmax = Inf;

        case "s_{div}"
            D("s_{div}") = makedist("Uniform",0.01,0.1);
            dmin = 0; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "r_{mig}"
            D("r_{mig}") = makedist("Uniform",1,3);
            dmin = 0; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "p_{lim}"
            D("p_{lim}") = makedist("Uniform",7.5*(1+eps()),16.5*(1-eps()));
            dmin = 0; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        otherwise
            error("have not defined a distribution for %s.",par_names(i))

    end

    D(par_names(i)) = truncate(D(par_names(i)),dmin,dmax);

end
