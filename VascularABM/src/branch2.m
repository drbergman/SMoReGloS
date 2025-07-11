%helper function to perform branching
%INPUT: current CapillaryList, and the segment where branching occurs
%OUTPUT: Updated CapillaryList with new capillary
%Version 5
%Date: 4/4/2012
%Last user: Will feb 3, 2012
%           Kerri Aug 1, 2013
%Need to add:  disallow collisions; add real findDirection for new branch
%segment, i.e branch method in segmentS
%Check SegmentMatrix before branching

function [CapillaryList, SegmentMatrix,voxelgrid] = branch2(CapillaryList,SegmentMatrix, inputseg,Agentmat,voxelgrid, XYZh, prolifval,pixelsize)
    %global prolifval
    %gridsize = size(voxelgrid.VEGF);
    %Check number of arguments
%     if inputseg.Node2(3) == 990
%         uhohdHDp = 42;
%     end
    if nargin < 7
        error('Too few args')
    end

%      disp('in branch2')
%      %disp('branched')
%                            disp(inputseg.Node1)
%                            disp(inputseg.Node2)
%                            disp('capillary- ')
%                            disp(num2str(inputseg.CapillaryNum))
        nonempty = CapillaryList(~cellfun('isempty',CapillaryList));
        capnum = length(nonempty) + 1;
        newCap = capillary;
        newSeg = segmentS2;
        newSeg.Node1 = inputseg.Node2;

        %define CycleLength N(24,2)/2
        hcycle =  prolifval + 2.*randn(1,1);
        newSeg.CycleLength = round(hcycle/2);
        
        %find nearest hypoxic cell instead 
        if isempty(XYZh)
            [dir_vec] = inputseg.findDirection(Agentmat, pixelsize);
        else
            hdists = pdist2(newSeg.Node1, 20*XYZh);
            [minval, minloc] = min(hdists);
            XYZmin = XYZh(minloc,:);
            dir_vec = (20*XYZmin - newSeg.Node1)/minval; %normalized vector
        end
        if any(dir_vec) %check if boundary
            tempnode = inputseg.Node2 + 3*ceil((dir_vec.*pixelsize));% ceil((dir_vec.*min_length));%V4 %We need this since we multiply vector by 4 %KERRI V2.3
%             if checkBoundaries(round(tempnode/pixelsize),gridsize)  %V4%WILL ADDED since *4 now, need to reduce new sprout length incrementally until not outside boundary
%                 for m = 1:min_length
%                     tempnode = inputseg.Node2 + ceil((dir_vec.*(min_length-m)));%V4
%                     if ~checkBoundaries(round(tempnode/pixelsize),gridsize)
%                         break;
%                     end
%                 end
%             end     

            temp = tempnode;

            %check that the new node is inside boundary
            if ~isempty(temp)
                gridsize = size(voxelgrid.Agent);
                [binarybound] = checkBoundaries3(round(temp/pixelsize), gridsize);
            end
            
            %Vessel regression
            SegmentMatrixdummy = SegmentMatrix;
            notmature = SegmentMatrix(:,9) == 0;
            SegmentMatrixdummy(notmature,:) = [];
%Put in MatureNodes instead
            %SegmentMatrixdummy = MatureNodes;
            binreg = 0;
            for mats = 1:size(SegmentMatrixdummy,1)
                
                %check Node2 and intersections
                %temp is the Node2
                if SegmentMatrixdummy(mats,7) == 1|| SegmentMatrixdummy(mats,7) == 2 || SegmentMatrixdummy(mats,7) == 3 || SegmentMatrixdummy(mats,7) == 4 || SegmentMatrixdummy(mats,7) == 5 || SegmentMatrixdummy(mats,7) == 6 || SegmentMatrixdummy(mats,7) == 7 || SegmentMatrixdummy(mats,7) == 8 
                    %don't include it
                else
                    r2 = SegmentMatrixdummy(mats,4:6);
                    r1 = SegmentMatrixdummy(mats,1:3);
                    v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector; = unit vector
                    %Define Node2
                    
                    
                    %ADDED KERRI V1.3
                    [axval axpos] =  min(1 - abs(v)); % find the axis closest to 1 to use
                    units = floor((r2(axpos)-r1(axpos))/v(axpos)); %The number of unit iterations
                    
                    if units > 0
                        for jj = 1:units
                            r3 = r1+jj.*v; %Define next position
                            %Is it within 100 of the Node2?
                                
                                if ((temp(1) - r3(1))^2 + (temp(2) - r3(2))^2 + (temp(3)-r3(3))^2) < 50^2
                                    binreg = 1;
                                end
                    
                        end %for jj
                    end %if untis
                    
                end %if not the first two capillaries
                
            end %for mats
            
           % binreg
            
            if ~isempty(temp) && binarybound == 0 && binreg == 0 %Only branch if not outside boundary and not regressed
              %  display('branching!');

                newSeg.Node2 = [(temp(1)),(temp(2)),(temp(3))];
                if pdist2(newSeg.Node1, newSeg.Node2,'euclidean') > 20
                    error('branch tip is too long')
                end
                newSeg.CapillaryNum = capnum; %set new segment's new capillary num
                newSeg.Listpos = 1; %set new segment at seg list position 1
                newSeg.updateDirection(); %update the direction vector
                %[voxelgrid.Agent] = fillAgentGrid(voxelgrid.Agent, newSeg.Node1, newSeg.Node2, newSeg.Radius);
                newCap.addSegment(newSeg); %add the segment to the capillary
                newSeg.activate; %activate the segment
                newSeg.Type = 1; %Define as tip cell
                hcycle =  prolifval + 2.*randn(1,1);
                newSeg.CycleLength = round(hcycle/2);
                
                 %fill in SegmentMatrix
                 epos = find(SegmentMatrix(:,1) == 0, 1 );
                  if isempty(epos)
                      epos = length(SegmentMatrix(:,1)) +1;
                  end
                 SegmentMatrix(epos,1:3) =  newSeg.Node1;
                 SegmentMatrix(epos,4:6) =  newSeg.Node2;
                 SegmentMatrix(epos,7) =  newSeg.CapillaryNum;
                 SegmentMatrix(epos,8) = newSeg.Listpos;
                 SegmentMatrix(epos,9) =  newSeg.Mature;
                
                %check whether the segment has to deactivate
%                 binarydeact = newSeg.canDeactivate(Agentmat, pixelsize);
%                 if binarydeact == 1
%                    % disp('deactivate tip')
%                     newSeg.deactivate; %activate the segment
%                 end
                
                input = inputseg.Listpos; %V3
                CapillaryList{capnum} = newCap; %add the capillary to the CapillaryList
                inputcap = inputseg.CapillaryNum;
                CapillaryList{inputcap}.BranchList(input) = 1; %Set the position of the branching seg to be "branched"
                %Do Salt and Pepper - NOT in capillaries
                if inputcap > 2
                    if input < length(nonempty)
                        CapillaryList{inputcap}.BranchList(input+1) = 1; 
                    end
                    if input > 1
                        CapillaryList{inputcap}.BranchList(input-1) = 1;
                    end
                end
                
                %Add the new seg to the agent grid
                [voxelgrid, newSeg] = fillAgentGridTortNodes(voxelgrid, newSeg, pixelsize);
                
            end %if ~is empty
            
        end
        
%     %Plot Branching
%          figure
%      hold on
%          p=patch(isosurface(voxelgrid.Agent==1,0));
%          set(p,'facecolor','cyan' ,'edgecolor', 'none');
%             daspect([1 1 1])
%             isonormals(voxelgrid.Agent==1,p)
%             view(3);
% 
%             camlight
% 
%             lighting gouraud      
%         title('Branch end')
    
end