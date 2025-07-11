%A function to do anastomosis - either tip or stalk
          function [voxelgrid,SegmentMatrix, newseg, caplist, activeList, MatureNodes] = doAnastSM(newseg, voxelgrid, SegmentMatrix, gridindex, caplist, activeList, pixelsize, MatureNodes)
          if nargin < 5 %Removed vector input V4.2 KERRI 
                 error('Too few args')
          end
          %display('Anastomosis')
         
              %check the cell index
              Cap = voxelgrid.CapInd(gridindex);
              Pos = voxelgrid.PosInd(gridindex);
              if Cap == 0
                  disp('doAnastSM - 0');
                  voxelgrid.Agent(gridindex) = 0;
                  return
              end
              %newseg
              segment = caplist{Cap}.SegmentList{Pos};
              %istip =  segment.Active;
              istip = isempty(caplist{Cap}.SegmentList{Pos+1});
              
              if istip == 1
                  %Do tip cell anastomosis
                  %Connect to Node2 of tip cell
                  %display('is tip')
                 
                  %Remove the old agentfill - Not NEEDED
                  %[voxelgrid] = emptyAgentGrid2(voxelgrid, newseg);
                  
                  %Assign newseg's new node
                  newseg.Node2 = segment.Node2;
                  
                  %Fill in new agentgrid
                  [voxelgrid, newseg] = fillAgentGridTortNodes(voxelgrid, newseg, pixelsize);
                  
                  %Turn off both segments
                  caplist{Cap}.SegmentList{Pos}.deactivate;
                  newseg.deactivate;
                  
                  %Prevent further activation
                  caplist{Cap}.SegmentList{Pos}.Branched = 1;
                  newseg.Branched = 1;
                  
                  %Make the two anastomosing capillaries mature
                  for cpos = 1:Pos
                      caplist{Cap}.SegmentList{cpos}.Mature = 1;
                      spos = find(SegmentMatrix(:,7) == Cap &  SegmentMatrix(:,8) == cpos);
                      SegmentMatrix(spos,9) = 1;
                      MatureNodes = cat(1,MatureNodes,caplist{Cap}.SegmentList{cpos}.Node1);
                      MatureNodes = cat(1,MatureNodes,caplist{Cap}.SegmentList{cpos}.Node2);
                  end
                  for cpos2 = 1:newseg.Listpos
                      caplist{newseg.CapillaryNum}.SegmentList{cpos2}.Mature = 1;
                      spos = find(SegmentMatrix(:,7) == newseg.CapillaryNum &  SegmentMatrix(:,8) == cpos2);
                      SegmentMatrix(spos,9) = 1;
                      MatureNodes = cat(1,MatureNodes,caplist{newseg.CapillaryNum}.SegmentList{cpos2}.Node1);
                      MatureNodes = cat(1,MatureNodes,caplist{newseg.CapillaryNum}.SegmentList{cpos2}.Node2);
                  end
                  
                   %error('is tip')
              elseif istip == 0
                  %Do stalk cell anastomosis
                  %Connect to place of overlap
                  %display('is stalk')
                  
                  %Remove the old agentfill - NOT NEEDED
                  %[voxelgrid] = emptyAgentGrid2(voxelgrid, newseg);
                  
                  %Assign newseg's new node
                  [subx suby subz] = ind2sub(size(voxelgrid.Agent), gridindex); 
                   newnode = cat(2, subx, suby, subz); 
                  newseg.Node2 = newnode.*pixelsize; %Added 8/1/13 Kerri
                  
                  %Fill in new agentgrid
                  [voxelgrid,newseg] = fillAgentGridTortNodes(voxelgrid, newseg, pixelsize);
                  
                  %Turn off both segments
                  caplist{Cap}.SegmentList{Pos}.deactivate;
                  newseg.deactivate;
                  
                  %Prevent further activation
                  caplist{Cap}.SegmentList{Pos}.Branched = 1;
                  newseg.Branched = 1;
                  %error('is stalk')
                  
                  %Make the two anastomosing capillaries mature
                  for cpos = 1:Pos
                      caplist{Cap}.SegmentList{cpos}.Mature = 1;
                      spos = find(SegmentMatrix(:,7) == Cap &  SegmentMatrix(:,8) == cpos);
                      SegmentMatrix(spos,9) = 1;
                      MatureNodes = cat(1,MatureNodes,caplist{Cap}.SegmentList{cpos}.Node1);
                      MatureNodes = cat(1,MatureNodes,caplist{Cap}.SegmentList{cpos}.Node2);
                  end
                  for cpos2 = 1:newseg.Listpos
                      caplist{newseg.CapillaryNum}.SegmentList{cpos2}.Mature = 1;
                      spos = find(SegmentMatrix(:,7) == newseg.CapillaryNum & SegmentMatrix(:,8) == cpos2);
                      SegmentMatrix(spos,9) = 1;
                      MatureNodes = cat(1,MatureNodes,caplist{newseg.CapillaryNum}.SegmentList{cpos2}.Node1);
                      MatureNodes = cat(1,MatureNodes,caplist{newseg.CapillaryNum}.SegmentList{cpos2}.Node2);
                  end
              else
                  error('istip has an invalid value')
              end
              
              if nargout < 4 %Removed vector input V4.2 KERRI 
                 error('Too few args out')
              end
              
          end %doAnast
