%removeSphere2
%This fills a discrete grid with sphere of radius R, around the point p1
%This file also stores the indices and segment type for each Agent
%This file changes floor and ceil to round
%Version 2.2
%Last Revision: Kerri April 25, 2012

function [grid] = removeSphere2(grid, p1, R) 
gs = size(grid.Agent);
for i = max(1,round(p1(1)) - R^2):min(gs(1),round(p1(1)) + R^2)
    for j = max(1,round(p1(2)) - R^2):min(gs(2),round(p1(2)) + R^2)
        for k = max(1,round(p1(3)) - R^2):min(gs(3),round(p1(3)) + R^2)
            if ((i - p1(1))^2 + (j - p1(2))^2 + (k-p1(3))^2) <= R^2
                grid.Agent(i,j,k) = false;
                %grid.Indices(i,j,k,1:2) = 0;
                grid.CapInd(i,j,k) = int16(0);
                grid.PosInd(i,j,k) = int16(0);
                %tempCell.Type = segtype;
                %grid.Cell(i,j,k) = tempCell;
            end
        end
    end
end

