function M = simPatient(M)

% this function simulates a single patient defined by the input M

M = finishParameterSetup_Patient(M);

%% initialize inputs
M.t = 0;
M.i = 1; % step index
M = initializeEnvironment(M);
M = initializeTracked(M);

%% initialize figure
if M.plot_pars.plot_fig
    M = initializeFigure(M);
end

%%
if M.save_pars.make_save
    M = finishParameterSetup_Saves(M);
    M = saveInitialModelData(M);
else
    M.save_pars.save_model_state = false;
    M.next_save_time = Inf;
end

if M.plot_pars.plot_fig && M.plot_pars.make_movie
    M = initializeMovie(M);
    warning("off",'MATLAB:audiovideo:VideoWriter:mp4FramePadded')
    writeVideo(M.vid,print(M.fig.handle,'-r100','-RGBImage'))
end

for ei = 1:M.events.n % event index
    M = finishParameterSetup_Event(M,M.events.times(ei)-M.t);
    M = simForward(M);

    extendTimeSeriesPlots(M,ei)

    switch M.events.event_index(ei)
        case Inf % end of simulation
            break;
        otherwise
            error("what even is this event??????")
    end
end

if M.save_pars.make_save
    M = saveFinalModelData(M);
end

if M.plot_pars.plot_fig && M.plot_pars.make_movie
    close(M.vid)
    warning("on",'MATLAB:audiovideo:VideoWriter:mp4FramePadded')
end
