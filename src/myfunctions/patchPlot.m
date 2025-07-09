function [p,l] = patchPlot(ax,X,Y,opts)

arguments
    ax
    X
    Y
    opts.Color = [0 0 0]
    opts.EdgeColor = 'none'
    opts.FaceAlpha = 0.2
    opts.LineStyle = "-"
    opts.LineWidth = 2
    opts.DisplayName = ''

    opts.min_val = -Inf
    opts.max_val = Inf

    opts.patchPlotCoordsOptions = []
end
% opts = defaultPatchPlotOptions;
% if nargin==4
%     opts = overrideDefaultOptions(opts,input_opts);
% end

ax.NextPlot = "add";
[x,y_mean,pc] = patchPlotCoords(X,Y,opts.patchPlotCoordsOptions);

if opts.min_val > -Inf
    pc{2} = max(opts.min_val,pc{2});
end

if opts.max_val < Inf
    pc{2} = min(opts.max_val,pc{2});
end

p = patch(ax,pc{1},pc{2},opts.Color,"EdgeColor",opts.EdgeColor,"FaceAlpha",opts.FaceAlpha,"DisplayName",opts.DisplayName);
l = plot(ax,x,y_mean,"Color",opts.Color,"LineStyle",opts.LineStyle,"LineWidth",opts.LineWidth,"DisplayName",opts.DisplayName);

end

% function default_options = defaultPatchPlotOptions
% 
% default_options.Color = [0 0 0];
% default_options.EdgeColor = 'none';
% default_options.FaceAlpha = 0.2;
% default_options.LineStyle = "-";
% default_options.LineWidth = 2;
% default_options.DisplayName = '';
% 
% default_options.min_val = -Inf;
% default_options.max_val = Inf;
% 
% default_options.patchPlotCoordsOptions = [];
% 
% end

