%checkBoundaries

function [binary] = checkBoundaries(point, gridsize)
%This function takes in a point (x,y,z) and a size of a grid and determines whether
%the point is outside the grid 
%Output: returns 0 if the point is inside the grid
%                1 if the point is outside the grid
%Input: point = the point to test
%       gridsize = the size of the grid.  We assume the grid starts at
%       position [1 1 1];

%Check number of arguments
if nargin < 2
       error('Too few args')
end

valid = [0 0 0];

if point(1) < 1 || point(1) > gridsize(1)
    valid(1) = 1;
end


if point(2) < 1 || point(2) > gridsize(2)
    valid(2) = 1;
end


if point(3) < 1 || point(3) > gridsize(3)
    valid(3) = 1;
end

if sum(valid) > 0
    binary = 1;
else
    binary = 0;
end
    