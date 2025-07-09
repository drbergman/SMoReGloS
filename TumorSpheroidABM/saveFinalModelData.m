function M = saveFinalModelData(M)

% saves final spatial data (if save_model_state) and saves the tracked
% field

if M.save_pars.save_model_state
    load(sprintf("data/sims/%s/output_%08d.mat",M.save_pars.sim_identifier,M.save_index-1),"time");
    if M.t>time % make new save
        M = saveModelData(M);
    end
end

tracked = M.tracked;

if ~isequal(M.save_pars.fields_to_keep,"all")
    fn = fieldnames(tracked);
    for i = 1:numel(fn)
        if ~any(M.save_pars.fields_to_keep==fn{i})
            tracked = rmfield(tracked,fn{i});
        end
    end
end

if M.save_pars.interpolate_tracked
    t_min = 1440*tracked.t;
    t = round(t_min);
    % make sure that all the time points really are on the minute (within tolerance)
    assert(all(abs(t_min-t)<1e-9))
    fn = fieldnames(tracked);
    for i = 1:length(fn)
        tracked.(fn{i}) = interp1(t,tracked.(fn{i}),M.save_pars.t_min);
    end
    tracked.t = M.save_pars.t_min/1440; % set this here to get rid of any rounding errors
end
save(sprintf("data/sims/%s/output_final",M.save_pars.sim_identifier),"tracked",'-v7.3')