function M = initializeMovie(M)

if ~isfield(M.save_pars,"sim_identifier")
    M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given

    while exist(sprintf("data/sims/%s",M.save_pars.sim_identifier),"dir") % just in case this directory already exists somehow (not sure how to processes could start at the same time to the millisecond and then one create this folder before the other looks for it)
        M.save_pars.sim_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
    end

    mkdir(sprintf("data/sims/%s",M.save_pars.sim_identifier))
end

if ~exist(sprintf("data/sims/%s/movie",M.save_pars.sim_identifier),"dir")
    mkdir(sprintf("data/sims/%s/movie",M.save_pars.sim_identifier))
end

M.vid = VideoWriter(sprintf("data/sims/%s/movie/out.mp4",M.save_pars.sim_identifier),'MPEG-4');
M.vid.FrameRate = 1/M.pars.max_dt;
open(M.vid)
