%checkAnastbranch
%Check whether the new segment hits anything and needs to anastomose
%This is specifically for segments that are the start of a branch

function [doanast, gridindex] = CheckAnastBranch(grid, newseg, capillarylist, pixelsize)
%This function compares the old grid to the newgrid with the filled
%segments.  It checks whether there is any overlap and if it overlaps with
%a segment other than the previous seg
%grid must have a grid.Agent and grid.Indices

% %MakeTempGrid
% tempgrid.Agent = false(size(grid.Agent));
% %tempgrid.Indices = zeros(size(grid.Indices));
% %Fill in Tempgrid
% [tempgrid, tempseg] = fillAgentGridTortNodes(tempgrid, newseg, pixelsize);
% %%Compare grids
% overlaplist = find(tempgrid.Agent & grid.Agent);
% % oldCap = oldseg.CapillaryNum;
% % oldPos = oldseg.Listpos;

[overlaplist,tempgrid] = makeOverlapList(newseg,pixelsize,grid);


%Make list of overlapping indices
CapList = grid.CapInd(overlaplist);
PosList = grid.PosInd(overlaplist);

next = 1;
OList = zeros(length(CapList),1);
for ii = 1:length(CapList)
    NodeOne = capillarylist{CapList(ii)}.SegmentList{PosList(ii)}.Node1;
    NodeTwo = capillarylist{CapList(ii)}.SegmentList{PosList(ii)}.Node2;
    %Are they the same nodes?
    N1Truth = NodeOne == newseg.Node1; 
    N2Truth = NodeTwo == newseg.Node1;
    if sum(N1Truth) < 3 && sum(N2Truth) <3 
        %Make list of actual overlapping ones
       % display('Anast Branch')
        OList(next) = ii;
        next = next +1;
    end
end

%Remove zeros
OList(OList == 0) = [];

%Find the ones on the list not equal to the oldseg pos/cap
%OList = find(PosList ~= oldPos & CapList ~= oldCap);

if ~isempty(OList)
    %display('Not empty')
    doanast = 1;
    %Find closest grid point to anastomose
    truelist = overlaplist(OList); 
    [subx suby subz] = ind2sub(size(tempgrid.Agent), truelist); 
    dataY = cat(2, subx, suby, subz); 
    dataX = round(newseg.Node1./pixelsize);
    Edist = pdist2(dataX, dataY);
    [minval minind] = min(Edist);
    gridindex = truelist(minind);
    %display('check anast branch')
else
    doanast = 0;
    gridindex = {};
   
end


