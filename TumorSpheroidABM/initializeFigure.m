function M = initializeFigure(M)

M.fig.handle = figure;
set(M.fig.handle,'units','normalized','outerposition',[0 0 1 1])

M.fig.nrows = 5;
M.fig.ncols = 5;

%% indices for the subplots in ax
M.fig.scatter_ind = 1;
M.fig.tum_density_ind = 2;
M.fig.population_ind = 3;
M.fig.event_ind = 4;
M.fig.subpop_ind = 5;
M.fig.tum_prob_ind = 6;
M.fig.cell_slice_ind = 7;

ind_names = fieldnames(M.fig);
ind_names = ind_names(endsWith(ind_names,'ind'));

n_axes = 0;
for i = 1:length(ind_names)
    n_axes = max(n_axes,M.fig.(ind_names{i}));
end
M.fig.ax = gobjects(n_axes,1);

%% set up axes
cellfun(@(f) evalin('caller',[f ' = M.fig.' f ';']), fieldnames(M.fig));

scatter_locs = [1,2*ncols+3];
cell_slice_locs = [3*ncols+1,4*ncols+2];

tumor_density_locs = 0*ncols + 4;

population_locs = 1*ncols + (4:5);
event_locs = 3*ncols + (4:5);
subpop_locs = 2*ncols + (4:5);

tum_probs_locs = 0*ncols + 4:5;


%% scatter plot
% phase_colors = [0.23,0.7,0.34;... % color for m phase
%                 0.05,0.06,0.58;... % color for g0
%                 0.95,0.98,0.06]; % color for g1
phase_names = ["G1","S","G2","M"];
phase_size = [30,40,60,60];
phase_colors = lines(4);
tumor_colors = winter(3);
tumor_colormap = flipud(winter);

if M.plot_pars.plot_location
    M.fig.ax(scatter_ind) = subplot(nrows,ncols,scatter_locs);
    M.fig.ax(scatter_ind).Box = 'off';
    M.fig.ax(scatter_ind).NextPlot = 'add';
    M.fig.ax(scatter_ind).Color = M.fig.ax(scatter_ind).Parent.Color;
    M.fig.ax(scatter_ind).XTick = [];
    M.fig.ax(scatter_ind).YTick = [];
    M.fig.ax(scatter_ind).XColor = 'none';
    M.fig.ax(scatter_ind).YColor = 'none';

    if M.setup.ndims==3
        M.fig.ax(scatter_ind).ZTick = [];
        M.fig.ax(scatter_ind).ZColor = 'none';

        for i = 1:4
            ind = M.tumor(:,M.I.phase)==i;
            M.fig.scatter_plots(i) = scatter3(M.tumor(ind,M.I.subs(1))-M.grid.center(1),M.tumor(ind,M.I.subs(2))-M.grid.center(2),M.tumor(ind,M.I.subs(3))-M.grid.center(3),phase_size(i),phase_colors(i,:),'o','filled',...
            'DisplayName',phase_names(i));
        end

        legend(M.fig.ax(M.fig.scatter_ind),'Location','SouthWest','AutoUpdate','off','Color',"none")

        M.fig.ax(M.fig.scatter_ind).Legend.Position(1) = M.fig.ax(M.fig.scatter_ind).Position(1) - M.fig.ax(M.fig.scatter_ind).Legend.Position(3);
        axis(M.fig.ax(M.fig.scatter_ind),[1,M.grid.size(1),1,M.grid.size(2),1,M.grid.size(3)] - repelem(M.grid.center,1,2))

        M.plot_pars.m = sqrt(sum(M.grid.center.^2));
        view(M.fig.ax(M.fig.scatter_ind), M.plot_pars.m*[0.8660,0,0.5])
        axis square
    else
        for i = 1:4
            ind = M.tumor(:,M.I.phase)==i;
            M.fig.scatter_plots(i) = scatter(M.tumor(ind,M.I.subs(1))-M.grid.center(1),M.tumor(ind,M.I.subs(2))-M.grid.center(2),phase_size(i),phase_colors(i,:),'o','filled',...
            'DisplayName',phase_names(i));
        end

        M.fig.ax(M.fig.scatter_ind).Position(1) = 0;
        legend(M.fig.ax(M.fig.scatter_ind),'Location','Northwest','AutoUpdate','off','Color',"white")

        axis(M.fig.ax(M.fig.scatter_ind),[1,M.grid.size(1),1,M.grid.size(2)] - repelem(M.grid.center,1,2))

        axis square

    end
end

%% slice plots
if M.setup.ndims==3
    warning("Do not color the slice for a 3d sim by phase yet")
    warning('off','MATLAB:contour:NonFiniteData')
    %% cell slice plot
    M.fig.ax(cell_slice_ind) = subplot(nrows,ncols,cell_slice_locs);
    M.fig.ax(cell_slice_ind).NextPlot = 'add';
    M.fig.ax(cell_slice_ind).XLim = [1,M.grid.size(1)];
    M.fig.ax(cell_slice_ind).YLim = [1,M.grid.size(2)];
    M.fig.cell_slice_plot(1) = scatter3([],[],[],60,tumor_colors(1,:),'o','filled',...
        'DisplayName','Tumor Cells','MarkerFaceAlpha',0.75);

    z_mid_log_tumor = M.tumor(:,M.I.subs(3))==M.grid.center(3);
    M.fig.cell_slice_plot(1).XData = M.tumor(z_mid_log_tumor,M.I.subs(1));
    M.fig.cell_slice_plot(1).YData = M.tumor(z_mid_log_tumor,M.I.subs(2));
    M.fig.cell_slice_plot(1).ZData = ones(length(M.fig.cell_slice_plot(1).XData),1);

    M.fig.ax(cell_slice_ind).Title.String = 'Cell Slice (z=z_{mid})';
    M.fig.ax(cell_slice_ind).XLabel.String = 'X coord (cells)';
    M.fig.ax(cell_slice_ind).YLabel.String = 'Y coord (cells)';

    axis square

    if ~M.plot_pars.plot_location
        legend(M.fig.cell_slice_plot,'Location','bestoutside','AutoUpdate','off')
    end

    %% density projection plots
    %% tumor density plot
    M.fig.ax(tum_density_ind) = subplot(nrows,ncols,tumor_density_locs);
    [~,M.fig.tum_density_plot] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
    M.fig.ax(tum_density_ind).Title.String = 'Tumor Projected onto (X,Y)';
    M.fig.ax(tum_density_ind).XLabel.String = 'X coord (cells)';
    M.fig.ax(tum_density_ind).YLabel.String = 'Y coord (cells)';
    M.fig.tum_density_colorbar = colorbar;

    if M.NT>1
        [N,~,~] = histcounts2(M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2)),(M.grid.xx(1)-.5):(M.grid.xx(end)+.5),(M.grid.yy(1)-.5):(M.grid.yy(end)+.5),'Normalization','countdensity');
        if all(size(N)>1)
            M.fig.tum_density_plot.ZData = N';
        else
            M.fig.tum_density_plot.ZData = [];
        end
    end

    colormap(M.fig.ax(tum_density_ind),tumor_colormap);

    axis square
end

%% time series plots
[xtick_vals,order] = sort(M.events.times);
symbols = {'I','^','.','C'};
vals = [0,1,2,Inf];
xx = [0 max(eps(),xtick_vals(1))];

[symbol_ind,~] = find((vals==M.events.event_index(order))');
xtick_labels = symbols(symbol_ind);

for i = length(xtick_vals):-1:2
    if xtick_vals(i)==xtick_vals(i-1)
        xtick_labels{i-1} = ';';
        xtick_labels(i) = [];
    end
end
xtick_vals = unique(xtick_vals);

%% population plot
M.fig.ax(population_ind) = subplot(nrows,ncols,population_locs);
M.fig.ax(population_ind).Title.String = 'Population Numbers';
M.fig.ax(population_ind).NextPlot = 'add';
M.fig.ax(population_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(population_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(population_ind).YLabel.String = 'Count';
M.fig.ax(population_ind).XLim = xx;
M.fig.ax(population_ind).Color = [0 0 0];

M.fig.population_plots(1) = plot(M.t,M.setup.N0,'w-','LineWidth',2,'DisplayName','Tumor');
% M.fig.population_plots(2) = plot(M.t,M.setup.NI0,'b-','LineWidth',2,'DisplayName','N_I');

legend(M.fig.ax(population_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% subpopulations proportion plot
M.fig.ax(subpop_ind) = subplot(nrows,ncols,subpop_locs);
M.fig.ax(subpop_ind).Title.String = 'Sub-Population Proportions';
M.fig.ax(subpop_ind).NextPlot = 'add';
M.fig.ax(subpop_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(subpop_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(subpop_ind).YLabel.String = 'Proportion';
M.fig.ax(subpop_ind).XLim = xx;
M.fig.ax(subpop_ind).YLim = [0 1];
M.fig.ax(subpop_ind).Color = [0 0 0];

for i = 1:M.cycle.n_phases
    M.fig.subpop_plots(i) = plot(M.t,M.tracked.phases(1,i)/M.NT,'Color',phase_colors(i,:),'LineStyle','-','LineWidth',2,'DisplayName',phase_names(i));
end

L = legend(M.fig.ax(subpop_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color);
L.Position(1) = sum(M.fig.ax(scatter_ind).Position([1,3]));

%% event plot
M.fig.ax(event_ind) = subplot(nrows,ncols,event_locs);
M.fig.ax(event_ind).Title.String = 'Event Rates in Update';
M.fig.ax(event_ind).XLabel.String = 'Days';
M.fig.ax(event_ind).YLabel.String = 'Rate (days^{-1})';
M.fig.ax(event_ind).NextPlot = 'add';
% M.fig.ax(event_ind).XAxis.TickValues = xtick_vals;
% M.fig.ax(event_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(event_ind).XLim = xx;
M.fig.ax(event_ind).YScale = 'log';
M.fig.ax(event_ind).Color = [0 0 0];

M.fig.events_plots(1) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Proliferation');
M.fig.events_plots(2) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','--','LineWidth',2,'DisplayName','Contact Inhibitions');
M.fig.events_plots(3) = plot(M.t,0,'Color',[0.9,0.1,0.1],'LineStyle','-','LineWidth',2,'DisplayName','Apoptosis');
M.fig.events_plots(4) = plot(M.t,0,'Color',[0.9,0.1,0.1],'LineStyle','--','LineWidth',2,'DisplayName','Chemo Death');

L = legend(M.fig.ax(event_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color);
L.Position(1) = sum(M.fig.ax(scatter_ind).Position([1,3]));

%% not time series
%% tumor probabilities
M.fig.ax(tum_prob_ind) = subplot(nrows,ncols,tum_probs_locs);
warning("Apoptosis for arrested cells may make this inappropriate. Not dealing with it now because we are not using this at the moment.")
max_tum_prob = sum(1-exp(-[max(M.cycle_pars.transition_rates),M.pars.apop_rate,M.pars.move_rate]*M.pars.max_dt));
M.fig.ax(tum_prob_ind).YLim = [0 max_tum_prob];
M.fig.ax(tum_prob_ind).Title.String = 'Tum Outcome Probabilities';
M.fig.ax(tum_prob_ind).NextPlot = 'replacechildren';
M.fig.ax(tum_prob_ind).XAxis.Visible = "off";

M.fig.tum_prob_bar = bar(NaN,NaN(3,1),1,'stacked','EdgeAlpha',0);

events = {'transition','apop','move'};
order = [3,2,1];
L = legend(M.fig.ax(tum_prob_ind),flip(M.fig.tum_prob_bar),flip(events(order)),'AutoUpdate','off','Location','northwest');
L.Position(1) = sum(M.fig.ax(scatter_ind).Position([1,3]));

%% set font size
for i = 1:length(M.fig.ax)
    if ishandle(M.fig.ax(i))
        M.fig.ax(i).FontSize = M.plot_pars.default_font_size;
    end
end

drawnow