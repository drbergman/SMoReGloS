%fillAgentGrid2
% This fills the  grid with ones in the approximate location of the agents,
%meaning the segment that extends from point r1 to r2 qith radius radius
%Uses MakeSphere2
%Version 3
%Last Revision: Kerri July 18, 2012

function [grid] = fillAgentGrid2(grid, seg, pixelsize)
%inputs: grid, the grid you will be filling with ones, r1 the first
%position of the segment, r2 the last position of the segment, radius the
%radius of the segment, indices = the capillary position and seg position
%of the segment

if nargin < 2 
    error('Too few args') 
end
%disp('In fill Agent')
%assign variables
r1 = round(seg.Node1./pixelsize); 
r2 = round(seg.Node2./pixelsize);
radius = round(seg.Radius/pixelsize);
indices = [seg.CapillaryNum seg.Listpos];

%Fill the first position of the segment
grid = makeSphere2(grid, r1, radius, indices); %or make Circle(grid, r3, radius, 3)

v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector; = unit vector
%ADDED KERRI V1.3
[axval axpos] =  min(1 - abs(v)); % find the axis closest to 1 to use
units = floor((r2(axpos)-r1(axpos))/v(axpos)); %The number of unit iterations

if units > 0   
    for j = 1:units 
        r3 = r1+j.*v; %Define next position
        grid = makeSphere2(grid, r3, radius, indices); %Set the sphere to 1
    end
end

%Fill the first position of the segment 
grid = makeSphere2(grid, r2, radius, indices);

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
