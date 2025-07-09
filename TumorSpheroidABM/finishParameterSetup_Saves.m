function M = finishParameterSetup_Saves(M)

% finished any necessary setup for save parameters, including making
% necessary output folders

if M.save_pars.dt < Inf
    M.save_pars.save_model_state = true;
    M.next_save_time = 0;
    M.save_index = 0;
else
    M.save_pars.save_model_state = false;
end

if ~isfield(M.save_pars,"sim_identifier")
    M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmm")); % default to this for determining an id if none given
    while exist(sprintf("data/sims/%s",M.save_pars.sim_identifier),"dir") % just in case this directory already exists somehow (not sure how to processes could start at the same time to the millisecond and then one create this folder before the other looks for it)
        M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmmss")); % default to this for determining an id if none given
    end
end

if isfield(M.save_pars,"idx_in_cohort") && ~endsWith(M.save_pars.sim_identifier,string(M.save_pars.idx_in_cohort))
    M.save_pars.sim_identifier = sprintf("%s_%d",M.save_pars.sim_identifier,M.save_pars.idx_in_cohort);
end

mkdir(sprintf("data/sims/%s",M.save_pars.sim_identifier))

max_grid_sub = max(M.grid.size);
max_grid_sub_log2 = ceil(log2(max_grid_sub+1));
if max_grid_sub_log2 <= 8
    M.save_pars.integrify = @uint8;
elseif max_grid_sub_log2 <= 16
    M.save_pars.integrify = @uint16;
elseif max_grid_sub_log2 <= 32
    M.save_pars.integrify = @uint32;
elseif max_grid_sub_log2 <= 64
    M.save_pars.integrify = @uint64;
else
    error("grid too large to store!")
end
