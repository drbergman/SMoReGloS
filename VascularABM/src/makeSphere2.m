%makeSphere2
%This fills a discrete grid with sphere of radius R, around the point p1
%This file also stores the indices and segment type for each Agent
%This file changes floor and ceil to round
%Version 2.3
%Last Revision: Daniel 2023

%% Daniel's update for Kerri's version using grid
function [grid] = makeSphere2(grid, p1, R, indices,onlyAgentGridNeeded)
if nargin==4
    onlyAgentGridNeeded = false;
end
if int16(indices(1)) == 0
    error('CapInd is 0!!!')
end
if isfield(grid,"sz")
    gs = grid.sz;
else
    gs = size(grid.Agent);
end
[X,Y,Z] = ndgrid(max(1,ceil(p1(1)-R)):min(gs(1),floor(p1(1)+R)),max(1,ceil(p1(2)-R)):min(gs(2),floor(p1(2)+R)),max(1,ceil(p1(3)-R)):min(gs(3),floor(p1(3)+R)));
D = (X-p1(1)).^2+(Y-p1(2)).^2+(Z-p1(3)).^2;
I = find(D<=R^2);
ind = sub2ind(gs,X(I),Y(I),Z(I));
grid.Agent(ind) = true;
if ~onlyAgentGridNeeded
    grid.CapInd(ind) = int16(indices(1));
    grid.PosInd(ind) = int16(indices(2));
end