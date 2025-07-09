function plotFunction_EndStep(M)

%% cell locations
if M.plot_pars.plot_location
    for i = 1:4
        ind = M.tumor(:,M.I.phase)==i;

        if M.setup.ndims == 3
            M.fig.scatter_plots(i).XData = M.tumor(ind,M.I.subs(1))-M.grid.center(1);
            M.fig.scatter_plots(i).YData = M.tumor(ind,M.I.subs(2))-M.grid.center(2);
            M.fig.scatter_plots(i).ZData = M.tumor(ind,M.I.subs(3))-M.grid.center(3);

            view(M.fig.ax(M.fig.scatter_ind), M.plot_pars.m*[cos(M.t*2*pi)*0.8660,sin(M.t*2*pi)*0.8660,0.5])
        else
            M.fig.scatter_plots(i).XData = M.tumor(ind,M.I.subs(1))-M.grid.center(1);
            M.fig.scatter_plots(i).YData = M.tumor(ind,M.I.subs(2))-M.grid.center(2);
        end
    end
end

if M.setup.ndims == 3
    %% slice
    z_mid_log_tumor = M.tumor(:,M.I.subs(3))==M.grid.center(3);
    M.fig.cell_slice_plot(1).XData = M.tumor(z_mid_log_tumor,M.I.subs(1));
    M.fig.cell_slice_plot(1).YData = M.tumor(z_mid_log_tumor,M.I.subs(2));
    M.fig.cell_slice_plot(1).ZData = ones(length(M.fig.cell_slice_plot(1).XData),1);

    view(M.fig.ax(M.fig.cell_slice_ind),2);

    %% projections
    if M.NT>1
        [N,~,~] = histcounts2(M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2)),(M.grid.xx(1)-.5):(M.grid.xx(end)+.5),(M.grid.yy(1)-.5):(M.grid.yy(end)+.5),'Normalization','countdensity');
        if all(size(N)>1)
            M.fig.tum_density_plot.ZData = N';
        else
            M.fig.tum_density_plot.ZData = [];
        end
    end

end

%% tracked quantities
%% population plot
M.fig.population_plots(1).XData(end+1) = M.t;
M.fig.population_plots(1).YData(end+1) = M.NT;

%% subpopulations proportion plot
for i = 1:M.cycle.n_phases
    M.fig.subpop_plots(i).XData(end+1) = M.t;
    M.fig.subpop_plots(i).YData(end+1) = M.tracked.phases(M.i,i)/(M.NT);
end

%% event plots
M.fig.events_plots(1).XData(end+1) = M.t;
M.fig.events_plots(1).YData(end+1) = M.tracked.tum_prolif(M.i) / (M.tracked.NT(M.i-1)*M.dt);
M.fig.events_plots(2).XData(end+1) = M.t;
M.fig.events_plots(2).YData(end+1) = M.tracked.tum_contact_inhibition(M.i) / (M.tracked.NT(M.i-1)*M.dt);
M.fig.events_plots(3).XData(end+1) = M.t;
M.fig.events_plots(3).YData(end+1) = M.tracked.tum_apop(M.i) / (M.tracked.NT(M.i-1)*M.dt);
M.fig.events_plots(4).XData(end+1) = M.t;
M.fig.events_plots(4).YData(end+1) = M.tracked.chemo_arrest(M.i) / (M.tracked.NT(M.i-1)*M.dt);

drawnow