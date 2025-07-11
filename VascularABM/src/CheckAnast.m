%checkAnast
%Check whether the new segment hits anything and needs to anastomose

function [doanast, gridindex] = CheckAnast(grid, newseg, oldseg, CapillaryList, pixelsize)
%This function compares the old grid to the newgrid with the filled
%segments.  It checks whether there is any overlap and if it overlaps with
%a segment other than the previous seg
%grid must have a grid.Agent and grid.Indices
if nargin < 5  
   error('Too few args')
end
            

% current_seed = rng;
% [overlaplist1,tempgrid1] = kerri_overlap_list(newseg,pixelsize,grid);
% rng(current_seed);
[overlaplist,tempgrid] = makeOverlapList(newseg,pixelsize,grid);

% assert(isequal(overlaplist1,overlaplist2))

% overlaplist = overlaplist1;

% zz = find(overlaplist==0);
% overlaplist(zz)=[];
oldCap = oldseg.CapillaryNum;
oldPos = oldseg.Listpos;

%Make list of overlapping indices
CapList = grid.CapInd(overlaplist);
PosList = grid.PosInd(overlaplist);

%Find the ones on the list not equal to the oldseg pos/cap
OList = find(PosList ~= oldPos | CapList ~= oldCap);
 %test that none of the other segment's nodes overlap with newcell Node1
 EiList = zeros(size(OList));
    for Ei = 1:length(OList)
        if CapList(OList(Ei)) == 0
            disp('CheckAnast - 0');
            %grid.Agent(OList(Ei)) = 0;
            continue
        end
        %check Node1
        NSN1 = round(newseg.Node1);
        ASN1 = round(CapillaryList{CapList(OList(Ei))}.SegmentList{PosList(OList(Ei))}.Node1);
        ASN2 = round(CapillaryList{CapList(OList(Ei))}.SegmentList{PosList(OList(Ei))}.Node2);
        check1 = sum(NSN1 == ASN1);
        check2 = sum(NSN1 == ASN2);
        if check1 == 3 || check2 == 3
            %The nodes are the same - Don't Anast
            EiList(Ei) = 1;
        else
            %Do Anast
            EiList(Ei) = 0;
        end
    end
%Remove the Node overlaps from OList    
EiBin = EiList == 1;
OList(EiBin) = [];

if ~isempty(OList)
    %display('Not empty')
    
    %Find closest grid point to anastomose
    truelist = overlaplist(OList); 
    [subx, suby, subz] = ind2sub(size(tempgrid.Agent), truelist); 
    dataY = cat(2, subx, suby, subz); 
    dataX = round(newseg.Node1/pixelsize);
    Edist = pdist2(dataX, dataY);
    %identify those less than 1 pixelsize away
    tooclose = Edist <= pixelsize;
    %remove from Edist & truelist
    Edist(tooclose) = [];
    truelist(tooclose) = [];
    
   
    if ~isempty(Edist)
        %error('Anast overlaps seg')
        [minval, minind] = min(Edist);
        doanast = 1;
    
        gridindex = truelist(minind);
        %newseg;
        %display('check anast')
    else
        doanast = 0;
        gridindex = {};
    end
    
    

    %error('same segs')
else
    doanast = 0;
    gridindex = {};
end

end

% function [overlaplist,tempgrid] = kerri_overlap_list(seg,pixelsize,grid)
% 
% %MakeTempGrid
% tempgrid.Agent = false(size(grid.Agent));
% %tempgrid.Indices = zeros(size(grid.Indices));
% %Fill in Tempgrid
% [tempgrid, tempseg] = fillAgentGridTortNodes(tempgrid, seg, pixelsize);
% %%Compare grids
% overlaplist = find(tempgrid.Agent & grid.Agent);
% end
