function cohort_pars = updateCohortTimer(cohort_pars,num_finished)

if mod(num_finished,cohort_pars.update_timer_every)==0
    v_n = toc(cohort_pars.batch_start);
    cohort_pars.mu_n = (((num_finished/cohort_pars.update_timer_every)-1)*cohort_pars.mu_n+2*v_n)/((num_finished/cohort_pars.update_timer_every)+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
    etr = cohort_pars.mu_n*(cohort_pars.n_batches-num_finished)/cohort_pars.update_timer_every;
    fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
        num_finished*cohort_pars.n_per_batch,cohort_pars.total_runs,100*num_finished*cohort_pars.n_per_batch/cohort_pars.total_runs,duration(0,0,v_n),duration(0,0,etr),duration(0,0,etr+toc(cohort_pars.start)))
    cohort_pars.batch_start = tic;
end
