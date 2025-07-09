function normalizeZLims(in)

% take a figure in and set all the axes to have the same z limits
% if input is set of axes, set them to all have the same z limits


switch class(in)
    case 'matlab.ui.Figure'
        h = findall(in,'type','axes');
        zls = arrayify(h,'ZLim');
        ZL = [min(zls(:,1)),max(zls(:,2))];
        set(h,'ZLim',ZL)

    case 'matlab.graphics.axis.Axes'
        zls = arrayify(in(:),'ZLim');
        ZL = [min(zls(:,1)),max(zls(:,2))];
        set(in,'ZLim',ZL)
end