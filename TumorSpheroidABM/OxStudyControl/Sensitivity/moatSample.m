function out = moatSample(x,M,par_names,D,nsamps,alpha,ci_relative_spread)

% runs nsamps at the point in parameter space defined by x. M has the base
% parameters of the simulation. par_names stores the correspondence between
% indices in x and parameters in M. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value.
for i = 1:numel(x)
    v_temp = icdf(D(par_names(i)),x(i));
    switch par_names(i)
        case "carrying_capacity"
            M.setup.carrying_capacity = v_temp;
        case "g1_to_s"
            M.cycle_pars.g1_to_s = v_temp;
        case "s_to_g2"
            M.cycle_pars.s_to_g2 = v_temp;
        case "g2_to_m"
            M.cycle_pars.g2_to_m = v_temp;
        case "m_to_g1"
            M.cycle_pars.m_to_g1 = v_temp;
        case "arrest_coeff_g1"
            M.chemo_pars.arrest_coeff_g1 = v_temp;
        case "arrest_coeff_g2"
            M.chemo_pars.arrest_coeff_g2 = v_temp;
        case "apop_rate"
            M.pars.apop_rate = v_temp;
        case "move_rate_microns"
            M.pars.move_rate_microns = v_temp;
        case "occmax_2d"
            M.pars.occmax_2d = min(7,floor(v_temp)); % make sure that this value cannot exceed 7 (otherwise cells will proliferate even when no space exists)

        otherwise
            error("Have not yet planned for %s to be varied.",par_names(i))
    end
end

M.save_pars.make_save = false; % make sure I'm not creating save files for any of this
M.plot_pars.plot_fig = false; % make sure I'm not plotting anything either

NT = zeros(nsamps,1);
for si = 1:nsamps
    temp = simPatient(M);
    NT(si) = temp.NT;
end

% std(NT)/sqrt(numel(NT)) gives the estimated SD of means of sample size numel(NT)
% tinv(1-0.5*alpha,numel(NT)-1) gives the t value for the 95% CI with numel(NT)-1 degrees of freedom
% Thus, the product of these (call it R) gives the radius of the CI about mean(NT) for estimating the true population mean, i.e. true ABM output on Day 3
% I want to insist that this radius, or spread, is less than 10% of the estimated mean, i.e. R/mean(NT) < 0.1, otherwise I continue taking more samples

relative_spread = tinv(1-0.5*alpha,numel(NT)-1)*std(NT)/(sqrt(numel(NT))*mean(NT));

if relative_spread > ci_relative_spread
    while relative_spread > ci_relative_spread % while the confidence interval for the mean is wider than 10% of the mean, keep going
        temp = simPatient(M);
        NT(end+1) = temp.NT;
        relative_spread = tinv(1-0.5*alpha,numel(NT)-1)*std(NT)/(sqrt(numel(NT))*mean(NT));
%         fprintf("  Spread after %d: %3.2f%%\n",numel(NT),100*relative_spread)
    end
else
end

out = mean(NT);