
function [overlaplist,tempgrid] = makeOverlapList(seg,pixelsize,grid)

radius = round(seg.Radius/pixelsize);
r1 = ceil(seg.Node1./pixelsize);
r2 = ceil(seg.Node2./pixelsize);
thedist = pdist2(r1,r2);
nz = ceil((thedist+2)*pi()*radius^2);
tempgrid.Agent = sparse([],[],[],numel(grid.Agent),1,nz);
tempgrid.sz = size(grid.Agent);
%tempgrid.Indices = zeros(size(grid.Indices));
%Fill in Tempgrid
[tempgrid, ~] = fillAgentGridTortNodes(tempgrid, seg, pixelsize, true);
% fprintf("Planned for %d. Got %d. Needed %d more nz.\n",nz,nnz(tempgrid.Agent),max(0,nnz(tempgrid.Agent)-nz))
%%Compare grids
overlaplist = find(tempgrid.Agent & grid.Agent(:));
