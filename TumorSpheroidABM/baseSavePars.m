function save_pars = baseSavePars()

save_pars.make_save = false; % whether or not to save any sim data
save_pars.dt = 1; % save by default once per day
save_pars.interpolate_tracked = false; % whether to thin out the tracked variables before saving at interpolated t values
save_pars.t_min = []; % time values (in minutes to help prevent rounding errors) at which to interpolate
save_pars.fields_to_keep = "all"; % which fields to keep; if set to "all", then keep all fields
save_pars.save_constants = true; % whether to save model constants by default
