function normalizeYLims(in)

% take a figure in and set all the axes to have the same y limits
% if input is set of axes, set them to all have the same y limits


switch class(in)
    case 'matlab.ui.Figure'
        h = findall(in,'type','axes');
        yls = arrayify(h,'YLim');
        YL = [min(yls(:,1)),max(yls(:,2))];
        set(h,'YLim',YL)

    case 'matlab.graphics.axis.Axes'
        yls = arrayify(in(:),'YLim');
        YL = [min(yls(:,1)),max(yls(:,2))];
        set(in,'YLim',YL)
end