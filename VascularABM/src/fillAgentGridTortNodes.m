%fillAgentGridToruosity
%This is modified from fillAgentGrid2
%This adds random sinusoidal toruosity to the network
% This fills the  grid with ones in the approximate location of the agents,
%meaning the segment that extends from point r1 to r2 qith radius radius
%Uses MakeSphere2
%Version 3
%Last Revision: Kerri August 20, 2015

function [grid, seg] = fillAgentGridTortNodes(grid, seg, pixelsize, onlyAgentGridNeeded)
%inputs: grid, the grid you will be filling with ones, r1 the first
%position of the segment, r2 the last position of the segment, radius the
%radius of the segment, indices = the capillary position and seg position
%of the segment

if nargin < 3
    error('Too few args') 
elseif nargin == 3
    onlyAgentGridNeeded = false;
end

%disp('In fill Agent')
%assign variables
r1 = ceil(seg.Node1./pixelsize);
r2 = ceil(seg.Node2./pixelsize);
radius = round(seg.Radius/pixelsize);
indices = [seg.CapillaryNum seg.Listpos];

%Find the nodes at length 5 away
thedist = pdist2(r1,r2);
steps = floor(thedist./10);
v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector; = unit vector
seg.NodeList = r1;
if steps > 1
    start = r1;
    for torts = 1:steps
        %go through the steps and add tort
    
        r3 = r1+ 10.*torts.*v; %Define next position
        r4 = abs(r3 +10*((rand(1,3)-0.5)));%Add random tortuosity - Dont want - numbers
        seg.updateNodes(r4);
        [grid] = fillAgentGrid3(grid, seg, start, r4, pixelsize, onlyAgentGridNeeded);
        start = r4;
    end
    %Fill the first position of the segment 
    [grid] = fillAgentGrid3(grid, seg, r4, r2, pixelsize, onlyAgentGridNeeded);
     seg.updateNodes(r2);
else
   %just go from p1 to p2
   [grid] = fillAgentGrid3(grid, seg, r1, r2, pixelsize,onlyAgentGridNeeded);
    seg.updateNodes(r2);
end


%Fill the first position of the segment 
%grid = makeSphere2(grid, r2, radius, indices);

% figure
%  p=patch(isosurface(grid==1,0));
%  %reducepatch(p, .5) 
% % nc = isocolors(slices==1,p);
%  set(p,'facecolor','cyan' ,'edgecolor', 'none');
%     daspect([1 1 1])
%     isonormals(grid==1,p)
%     view(3);
%     axis([20 50 10 80 40 400])
% 
%     camlight
% 
%     lighting gouraud 
