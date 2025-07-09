function normalizeXLims(f)

% take in a figure f and set all the axes to have the same x limits

h = findall(f,'type','axes');
xls = arrayify(h,'XLim');
XL = [min(xls(:,1)),max(xls(:,2))];
set(h,'XLim',XL)