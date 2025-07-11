classdef segmentS2 < handle
    %SEGMENTS = SegmentSimple This is a piece of the Capillary (to do)
    %   Segments are what the endothelial cells and capillaries consist of.
    %   Inheritance: They consist of different types: tip cells, and
    %   quiescent cells.  Default: quiescent. = Type
    %   Properties: They have a radius, length, and two nodes (1&2). They can 
            %be active, default (inactive), tip = the adjacent tip cell
            %(default = none), and D114 status (on or off) default = off. 
    %   Methods: They can become activated, proliferate,and
    %   branch depending on what type they are.
    %   D114 now affects: 1) number of activated segments
    %                     2) proliferation probability
    %   Includes a cell clock
    %   Includes tip cell, stalk cell and pushing stalk cell proliferation
    %   This uses tumor grid instead of voxelgrid
   
    
    %Last version = segmentS2; 
    %Version 7
    %Last revision: Will - FEB 10, 2012
    %               Kerri - May 16, 2012 
    %               Kerri - Oct 21 2013
    
     
    properties
        Type = 0; %Whether it is 0:quiescent (default), 1:tip cell, 2:start of endothelial cell 
        Radius = 2; %the radius (in microns) of the segment, default = 2, tip = 1;
        Length = 1; %length (in microns) of the segment, default = 0, tip = 5;
        Node1 = [1,1,1]; %position (x,y,z) of first node
        Node2 = [2,2,2]; %position (x,y,z) of second node %The leading node!
        Active = 0; %Whether the cell is active(a tip cell) or not - binary, default = 0 = inactive  
        Mature = 0; %Whether the segment is mature - due to Anast or hitting end (Mature segs cant prolif)
        D114 = 0; %Whether it is a D114 mutant 0 = no , 1 = yes
        Volume = 0; %Default stalk segment volume
        PorE = 0; %Whether the cell proliferated (0) or elongated (1) last iter (deafult=0) 
        lp = 0; %Either lpstalk or lptip = length changes due to prolif
        Listpos = 0; %The Cell's position in the segment container [index in current capillary, which capillary]
        direction = [1 1 1]; %This is for findDirection
        CapillaryNum = 0; % keeps track of which capillary this segment is in
        Branched = 0; %Whether this segment has branched before if so = quiescent
        LastAction = 0; %Last action of this segment 0 = prolif, 1 = migrate; ADDED WILL 2/10
        Pepper = 0; %Salt and pepper--pepper = inhibited, forced to stalk phenotype
        CycleLength = 12; %The length of thr cell cycle in iterations (~12 its = 24 hours)
        CellClock = 0; %The place in the cell cycle (in its) 
        EndoLength = 1; %Length of the endo cell up until this segment
        NodeList;
    end
    
    methods
        
         function [obj] = updateDirection(obj)
             obj.direction = (obj.Node2 - obj.Node1)./norm(obj.Node2 - obj.Node1);
         end
         function [obj] = updateNodes(obj,newnode)
             obj.NodeList = cat(1,obj.NodeList,newnode);
         end
 
         function [obj] = activate(obj)
%             if nargin < 2 
%                 error('Too few args') 
%             end
            obj.Active = 1; %Turn to active;
           % disp('Activated Segment');
%             activeList(obj.Listpos) = 1;
            
%             if nargout == 0
%                 error('Need to output activelist')
%             end
            
         end
         
         
         function [obj] = deactivate(obj)
%             if nargin < 2 
%                 error('Too few args') 
%             end
            obj.Active = 0; %Turn to deactive;
%             activeList(obj.Listpos) = 0;
%             
%             if nargout == 0
%                 error('Need to output activelist')
%             end
            
         end
         
         function [binary] = canDeactivate(obj,Agentmat, pixelsize)
             %global vegf_constant;
             if nargin < 2 
                error('Too few args') 
             end

             %Check whether the VEGF conditions are satisfied
             %Now check an area +- 1 for a cell
             %binary = 0;
             N1position = ceil(obj.Node1./20); %position of Node1 BSRI Remove pixelsize
             N2position = ceil(obj.Node2./20); %position of Node2 BSRI Remove pixelsize
             %grid
             VEGF1 = vegf_constant*sum(sum(Agentmat(max(N1position(1)-1,1):min(N1position(1)+1, size(Agentmat,1)),max(N1position(2)-1,1):min(N1position(2)+1, size(Agentmat,2)),max(N1position(3)-1,1):min(N1position(3)+1, size(Agentmat,3))))); %VEGF value at Node1
             VEGF2 = vegf_constant*sum(sum(Agentmat(max(N2position(1)-1,1):min(N2position(1)+1, size(Agentmat,1)),max(N2position(2)-1,1):min(N2position(2)+1, size(Agentmat,2)),max(N2position(3)-1,1):min(N2position(3)+1, size(Agentmat,3))))); %VEGF value of Node2
             if (VEGF1 < .05) || (VEGF2 < .05) %One node must have VEGF < .05 ng/ml (10% activation)
                 binary = 1;
             else
                 binary = 0;
             end
         end

         
         function [binary] = canActivate(obj, Agentmat, pixelsize)
              global vegf_constant;
              if nargin < 3 
                error('Too few args') 
              end
              
              %only do if segment has not branched before
              if obj.Branched == 1 %ADDED V3.3 with Will
                  binary = 0;
              else
                  %voxelgrid is the species grid
                  %numactive is the number of active tips is the capillary

                  %find the maximum number of active segments in each capillary
                  %depends on the cell's D114 status
%                   switch obj.D114
%                       case 0
%                           maxactive = 1;
%                       case 1
%                           maxactive = 2;
%                   end
% 
%                   %Check whether the VEGF conditions are satisfied
%                   numactive = sum(CapillaryList{obj.CapillaryNum}.ActiveList);
%                   if numactive < maxactive
                     %binary = 0;
                     N1position = ceil(obj.Node1./20); %position of Node1 BSRI Remove pixelsize
                     N2position = ceil(obj.Node2./20); %position of Node2 BSRI Remove pixelsize
                     %grid
                     VEGF1 = vegf_constant*sum(sum(sum(Agentmat(max(N1position(1)-1,1):min(N1position(1)+1, size(Agentmat,1)),max(N1position(2)-1,1):min(N1position(2)+1, size(Agentmat,2)),max(N1position(3)-1,1):min(N1position(3)+1, size(Agentmat,3)))))); %VEGF value at Node1
                     VEGF2 = vegf_constant*sum(sum(sum(Agentmat(max(N2position(1)-1,1):min(N2position(1)+1, size(Agentmat,1)),max(N2position(2)-1,1):min(N2position(2)+1, size(Agentmat,2)),max(N2position(3)-1,1):min(N2position(3)+1, size(Agentmat,3)))))); %VEGF value of Node2
           
                     if (VEGF1 > .5) && (VEGF2 > .5) %Both nodes must have VEGF > .5 ng/ml
                         binary = 1;
                     else 
                         binary = 0;
                     end
%                   else
%                      binary = 0;
%                   end
              end %only do if segment has not branched before
         end
         
         function [binaryP] = canProliferate(obj)
             %This is for tip cells only
             if pdist2(obj.Node2, obj.Node1) > 2 
                        binaryP = 1;
                % end
             else
                 binaryP = 0;
             end
                 
         end
    
           function [binaryP] =  canProliferateTip(obj)
            if obj.Active == 1 && pdist2(obj.Node2, obj.Node1) > 2  %Cell must be active (not quiescent)
                
              %depends on the cell's D114 status
              switch obj.D114
                  case 0 %Normal
                      %disp( 'case 0')
                      if rand < 1 %.03 ALWAYS PROLIFERATES!!!!!!!!
                          binaryP = 1;
                      else
                          binaryP = 0;
                      end
                  case 1
                      %disp('case 1')
                      if rand < 1 %.08
                          binaryP = 1;
                      else
                          binaryP = 0;
                      end
              end
            else
                binaryP = 0;
            end
             
         end
         
            
          %Second Stalk cell proliferation  %Aug 7
         function [newSeg] = proliferateStalk(obj, Listposval, prolifval)
             if nargin < 3 %Removed vector input  
                 error('Too few args')
             end
             
             %currently we will make the new cell the same size and
             %direction as previous cell
             %vector = obj.direction;
             vector2 = obj.Node2-obj.Node1;
             if sum(vector2 == 0) == 3
                error('proliferateStalk vector is 0')
            end
             
            %Make new segment
            newSeg = segmentS2;
            newSeg.Listpos = Listposval;
            
            %Define new Nodes
            newSeg.Node1 = obj.Node2; 
            newSeg.Node2 = round(obj.Node2+vector2); %New node %Will had vector*4 instead 
            
            %Define its capillaryNum %KERRI ADDED V2.2
            newSeg.CapillaryNum = obj.CapillaryNum; 
            
            %Change LastAction %Will ADDED V 3.3
            newSeg.LastAction = 0;
            
            %define CycleLength N(24,2)/2
            
            hcycle =  prolifval + 2.*randn(1,1);
            newSeg.CycleLength = round(hcycle/2);
            obj.CellClock = 0; %reset cell clock
            
            %Define new direction %KERRI ADDED V6.2 7/31/13
            newSeg.updateDirection;
            
            if nargout < 1
                error('Need to output newcell')
            end
            
         end
         
         %Tip and Stalk Cell Proliferation %V
         %This results in the addition of a tip segment 
         function [newSeg, newStalk, newEndo] = proliferate(obj, prolifval, pixelsize)
             %Note: input tip cell instead of stalk in here
             if nargin < 3 %Removed vector input V4.2 KERRI 
                 error('Too few args')
             end
             
             %global prolifval
             celltype = obj.Type;
             Listposval = obj.Listpos;
             %Listpos1 = obj.Listpos;
             %value for proliferation cell cycle (in hours)
             %prolifval = 42; - NOW IN CONTROL
            
            %Check whether it is the sprout or not
             if Listposval == 1 
                 %Do tip cell proliferation
               %  display('Doing tip cell prolif')
                 %define direction - Removed V4.2 KERRI
                 %VDir = vector;
                 vector2 = pixelsize*obj.direction; %Find the unit vector of direction
                 if sum(vector2 == 0) == 3
                     error('proliferate vector is 0')
                 end
                 
                 %Make new segment
                 newSeg = segmentS2;
                 newSeg.Listpos = Listposval + 1;
                 
                 %Define new Nodes
                 newSeg.Node1 = obj.Node2;
                 newSeg.Node2 = round(obj.Node2+vector2); %New node %Will had vector*4 instead
                 
                 if newSeg.Node1 == newSeg.Node2
                     error('prolif Nodes are the same')
                 end
                 
                 %Define EndoLength for new tip and stalk
                 % newSeg.EndoLength = pdist2(obj.Node1,obj.Node2); %tip
                 newSeg.EndoLength = simpleDist(obj.Node1,obj.Node2); %tip
                 % assert(isequal(temp,newSeg.EndoLength))
                 %obj.EndoLength = pdist2(newSeg.Node1,newSeg.Node2); %stalk
                 
                 %Make cell active and tip
                 newSeg.activate;
                 newSeg.Type = 1;
                 newSeg.CellClock = 0;
                 
                 %Deactivate old cell
                 newStalk = obj;
                 newStalk.deactivate;
                 newStalk.Type = 0;
                 newStalk.CellClock = 0;
                 
                 %Define its capillaryNum %KERRI ADDED V2.2
                 newSeg.CapillaryNum = obj.CapillaryNum;
                 
                 %Change LastAction %Will ADDED V 3.3
                 %newSeg.LastAction = 0;
                 
                 %define CycleLength N(24,2)/2
                 %prolifval
                 hcycle =  prolifval + 2.*randn(1,1);
                 newSeg.CycleLength = round(hcycle/2);
                 
                 %Define new direction %KERRI ADDED V4 - changed V5
                 newSeg.direction = obj.direction;
                 
                 %newStalk = {};
                 newEndo = 0;
             else
                 %Do stalk cell proliferation %March 12, 2013
               %  display('Doing stalk cell prolif')
                 %obj = tip cell
                 % vector = obj.Node2-obj.Node1; %Find the unit vector of direction 
                 vector2 = pixelsize*obj.direction;
                 if sum(vector2 == 0) == 3
                     error('proliferate vector is 0')
                 end
                 
                 %Make new stalk cell segment
                 newStalk = segmentS2;
                 newStalk.Listpos = Listposval; %stalk cell
                 newStalk.Type = 0;
                 %Define tip cell
                 newSeg = obj;
                 newSeg.Listpos = Listposval+1; %tip cell 
                 newSeg.Type = 1;
                 
                 %Check if you have a new endo cell
                 EndoSoFar = obj.EndoLength;
                 newStalk.EndoLength = EndoSoFar; %Update StalkCell EndoLength
                 
                     %Define new Nodes -stalk same as old tip
                 newStalk.Node1 = obj.Node1;
                 newStalk.Node2 = obj.Node2; %New node %Will had vector*4 instead
                 newStalk.NodeList = obj.NodeList;
                
                     % segL = pdist2(newStalk.Node1, newStalk.Node2);
                     segL = simpleDist(newStalk.Node1, newStalk.Node2);
                     % assert(isequal(segL,segLtemp))
                     EndoTotal = EndoSoFar + segL;
                     
                     
                     if EndoTotal > 30
                         newStalk.Type = 2;
                         newEndo = 1; %to record new cell location
                         %Reset EndoLength
                         newStalk.EndoLength = 0; %Update StalkCell EndoLength
                         newSeg.EndoLength = segL;
                     else
                         newSeg.EndoLength = EndoTotal; %Update Tip Cell EndoLength
                         newEndo =0;
                     end
                     
                      %Define new Nodes - tip is unit vector in direction of old tip
                 newSeg.Node1 = obj.Node2;
                 newSeg.Node2 = round(obj.Node2+vector2); %New node %Will had vector*4 instead
                 newSeg.NodeList = [];
                 
                 if newStalk.Node1 == newStalk.Node2
                     error('stalk prolif Nodes are the same')
                 end
                 
                 if newSeg.Node1 == newSeg.Node2
                     error('tip prolif Nodes are the same')
                 end
                 
                 %Make cell active
                 newSeg.activate;
                 
                 %Deactivate old cell
                 newStalk.deactivate;
                 
                 %Define its capillaryNum %KERRI ADDED V2.2
                 newStalk.CapillaryNum = obj.CapillaryNum;
                 newSeg.CapillaryNum = obj.CapillaryNum;
                 
                 %Update the new cycle length
                 hcycle =  prolifval + 2.*randn(1,1);
                 newStalk.CycleLength = round(hcycle/2);
                 newStalk.CellClock = 0;
                 newSeg.CellClock = 0;
                 
                 %Define new direction %KERRI ADDED V4 - changed V5
                 newStalk.updateDirection;
                 newSeg.updateDirection;
                 
                 %Change obj into stalkcell
                 %obj = newStalk;
             end %End of tip or stalk cell proliferation
             
            if nargout < 1
                error('Need to output newcell')
            end
            
         end
         
         %Can I migrate? NEED TO ADD OTHER CONDITIONS...NOT USED FOR NOW
         %Returns binary 1=yes, 0=no
         function [binaryM] = canMigrate(obj)
            if obj.LastAction < 4
                binaryM = 1;
            else
                binaryM = 0;
            end
         end
         
         %Method to perform migration. Is this even necessary...
         %This results in the migrating cell elongating towards high vegf
         function [obj] = migrate(obj, vector, elongation)
            %display('Migrating')
            if isnan(obj.Node1 + round(vector*elongation)) %Error check V5
                error('New node is NAN')
            end
            
%             if elongation > 30 
%                 error('elongation is > 30')
%             end
            
            obj.Node2 = obj.Node1 + round(vector*elongation);
            obj.updateDirection();
            
%             if pdist2(obj.Node1, obj.Node2) > 30 
%                 error('migration Nodes is > 30')
%             end
            
            if sum(obj.Node2 == obj.Node1) == 3
                error('Nodes are the same')
            end
           % obj.LastAction = obj.LastAction + 1;
         end
         
         %find direction based on nearest hypoxic cell
        function [dir_vec]= findDirectionH(obj,Agentmat,XYZh, pixelsize)
           % nonempty = CapillaryList(~cellfun('isempty',CapillaryList));
           % capnum = length(nonempty) + 1;
           % newCap = capillary;
           % newSeg = segmentS2;
           % newSeg.Node1 = inputseg.Node2;
           binh = 0;
            
            %find nearest hypoxic cell instead
            if isempty(XYZh)
                [dir_vec] = obj.findDirection(Agentmat, pixelsize);
            else
                % hdists = pdist2(obj.Node1, 20*XYZh);
                hdists = simpleDist(obj.Node1, 20*XYZh);
                % assert(isequal(hdiststemp,hdists))
                [minval, minloc] = min(hdists);
                XYZmin = XYZh(minloc,:);
                dir_vec = (20*XYZmin - obj.Node1)/minval; %normalized vector
                if (dir_vec(1) == 0 && dir_vec(2) == 0 && dir_vec(3) == 0) ||  dot(dir_vec, obj.direction) < 0%set self to neg infinity
                                binh = 1;
                end
            end
            
            %don't allow the tip to move backward
            if binh == 1
                [dir_vec] = obj.findDirection(Agentmat, pixelsize);
            end
        end
        
         %This is for unrestricted search
         function [dir_vec]= findDirection(obj,Agentmat, pixelsize) %to find largest VEGF gradient for migration. Returns normalized direction vec, target coordinates
            global vegf_constant;
            tipNode=obj.Node2;
            x = ceil(obj.Node2(1)/20);
            y = ceil(obj.Node2(2)/20);
            z = ceil(obj.Node2(3)/20);
            
            vegfbin = 3*Agentmat(x,y,z) ~= 0;
                            
            tipVEGF = vegf_constant*vegfbin; %VEGF at tip location
            
            coordinates = cell(1,28,3); % to store coordinates, vegf diff, and dir_vecs of each searched voxel for easy access
            
            ij=1;
            
            %Need to check that the startpoint and startpoint + 2 are not
            %outside of the grid
            startPoint = [x-1,y-1,z-1];
            thePoint = [x,y,z];
            gridsize = size(Agentmat);
            binary1 = checkBoundaries(startPoint, gridsize);
            binary2 = checkBoundaries(startPoint+2, gridsize);
            %check =  meshgrid([x-1,x,x+1], [y-1,y,y+1], [z-1,z,z+1]);
            %dir_vec = [];
            
            if binary1 == 0 && binary2 == 0  %ADDED KERRI V2.2
               [AAx,ABx,BBx] =  meshgrid([-1,0,1], [-1,0,1], [-1,0,1]);
                for ij = 1:numel(AAx)
                    coord = [AAx(ij), ABx(ij),BBx(ij)];
                            tempCoord = coord;
                            
                            vegfbin = Agentmat(x+tempCoord(1),y+tempCoord(2),z+tempCoord(3)) ~= 0;
                            tempVEGF = vegf_constant*vegfbin;
                            
                            dir_vec = tempCoord;

                            if (dir_vec(1) == 0 && dir_vec(2) == 0 && dir_vec(3) == 0) %set self to neg infinity
                                coordinates{1,ij,2}=-Inf;
                                continue;
                            end
                             if dot(dir_vec, obj.direction) < 0 %set points "behind" to neg infinity
                                coordinates{1,ij,2}=-Inf;
                                continue;
                            end
          
                            coordinates{1,ij,2} = tempVEGF - tipVEGF; 
                            coordinates{1,ij,1} = [x,y,z];
                            coordinates{1,ij,3} = dir_vec;
                              
                end

                %Set self-comparison to negative infinity
            
                maxindices = find([coordinates{1,:,2}] == max([coordinates{1,:,2}]));
                if length(maxindices) == 1
                    new_coord = coordinates{1,maxindices(1),1};
                    dir_vec= coordinates{1,maxindices(1),3};
%                     n_vec = dir_vec - thePoint;
%                     new_vec = tipNode+n_vec;

                else %choose randomly from maxes, if more than one. May need to add persistence bias!!
                    %display('random!');
                    random = randi(length(maxindices));
                    new_coord = coordinates{1,maxindices(random),1};
                    dir_vec= coordinates{1,maxindices(random),3};
%                      n_vec = dir_vec - thePoint;
%                     new_vec = tipNode+n_vec;

                end
            elseif binary1 == 0
                [AAx,ABx,BBx] =  meshgrid([-1,0], [-1,0], [-1,0]);
                for ij = 1:numel(AAx)
                    coord = [AAx(ij), ABx(ij),BBx(ij)];
                            tempCoord = coord;
                            
                            vegfbin = Agentmat(x+tempCoord(1),y+tempCoord(2),z+tempCoord(3)) ~= 0;
                            tempVEGF = vegf_constant*vegfbin;
                            
                            dir_vec = tempCoord;

                            if (dir_vec(1) == 0 && dir_vec(2) == 0 && dir_vec(3) == 0) %set self to neg infinity
                                coordinates{1,ij,2}=-Inf;
                                continue;
                            end
                            
                             if dot(dir_vec, obj.direction) < 0 %set points "behind" to neg infinity
                                coordinates{1,ij,2}=-Inf;
                                continue;
                            end
          
                            coordinates{1,ij,2} = tempVEGF - tipVEGF; 
                            coordinates{1,ij,1} = [x,y,z];
                            coordinates{1,ij,3} = dir_vec;
                              
                end

                %Set self-comparison to negative infinity
            
                maxindices = find([coordinates{1,:,2}] == max([coordinates{1,:,2}]));
                if length(maxindices) == 1
                    new_coord = coordinates{1,maxindices(1),1};
                    dir_vec= coordinates{1,maxindices(1),3};
%                      n_vec = dir_vec - thePoint;
%                     new_vec = tipNode+n_vec;

                else %choose randomly from maxes, if more than one. May need to add persistence bias!!
                    %display('random!');
                    random = randi(length(maxindices));
                    new_coord = coordinates{1,maxindices(random),1};
                    dir_vec= coordinates{1,maxindices(random),3};
%                      n_vec = dir_vec - thePoint;
%                     new_vec = tipNode+n_vec;

                end
            elseif binary2 == 0    
                 [AAx,ABx,BBx] =  meshgrid([0, 1], [0,1], [0,1]);
                for ij = 1:numel(AAx)
                    coord = [AAx(ij), ABx(ij),BBx(ij)];
                            tempCoord = coord;
                            
                            vegfbin = Agentmat(x+tempCoord(1),y+tempCoord(2),z+tempCoord(3)) ~= 0;
                            tempVEGF = vegf_constant*vegfbin;
                            
                            dir_vec = tempCoord;

                            if (dir_vec(1) == 0 && dir_vec(2) == 0 && dir_vec(3) == 0) %set self to neg infinity
                                coordinates{1,ij,2}=-Inf;
                                continue;
                            end
                            if dot(dir_vec, obj.direction) < 0 %set points "behind" to neg infinity
                                coordinates{1,ij,2}=-Inf;
                                continue;
                            end
          
                            coordinates{1,ij,2} = tempVEGF - tipVEGF; 
                            coordinates{1,ij,1} = [x,y,z];
                            coordinates{1,ij,3} = dir_vec;
                              
                end

                %Set self-comparison to negative infinity
            
                maxindices = find([coordinates{1,:,2}] == max([coordinates{1,:,2}]));
                if length(maxindices) == 1
                    new_coord = coordinates{1,maxindices(1),1};
                    dir_vec= coordinates{1,maxindices(1),3};
%                      n_vec = dir_vec - thePoint;
%                     new_vec = tipNode+n_vec;

                else %choose randomly from maxes, if more than one. May need to add persistence bias!!
                    %display('random!');
                    random = randi(length(maxindices));
                    new_coord = coordinates{1,maxindices(random),1};
                    dir_vec= coordinates{1,maxindices(random),3};
%                      n_vec = dir_vec - thePoint;
%                     new_vec = tipNode+n_vec;

                end
            else
                %ADDED KERRI V2.2
                dir_vec = []; %There is no vector - search could not be completed 
                
            end
            
         end
          
          function [dir_vec]= findDirectionPM(obj,voxelgrid,sigma,sigma_VEGF,gamma,sigma_g,C_gradient, pixelsize) %to find largest VEGF gradient for migration. Returns normalized direction vec, target coordinates
            tipNode=obj.Node2;
            x = round(obj.Node2(1)/pixelsize);
            y = round(obj.Node2(2)/pixelsize);
            z = round(obj.Node2(3)/pixelsize);
            tipVEGF = voxelgrid.VEGF(x,y,z); %VEGF at tip location
            
            coordinates = cell(1,98,3); % to store coordinates, vegf diff, and dir_vecs of each searched voxel for easy access
            
            %i=1;
%             sigma = pi/6;%(pi/4)*(tipVEGF.^3)./(.53^3+tipVEGF.^3); %Sensitivity to persistence (smaller is more persistent)
%             sigma_VEGF = .01; %Sensitivity of local VEGF search (smaller is small st. dev)
%             gamma = 1;
%             sigma_g = pi/12;
%             C_gradient = 1;
            
            %Need to check that the startpoint and startpoint + 2 are not
            %outside of the grid
            startPoint = [x-2,y-2,z-2];
            gridsize = size(voxelgrid.VEGF);
            binary1 = checkBoundaries(startPoint, gridsize);
            binary2 = checkBoundaries(startPoint+4, gridsize);
            
            if binary1 == 0 && binary2 == 0  %ADDED KERRI V2.2
                for i = 1:length(voxelgrid.Check)
                            tempCoord = [voxelgrid.Check(i,:)];
                            
                            tempVEGF = voxelgrid.VEGF(x+tempCoord(1), y+tempCoord(2), z+tempCoord(3));
                            
                            dir_vec = tempCoord;
                            persistence = (1/exp(1/sigma^2))*exp(dot(dir_vec./norm(dir_vec),obj.direction./norm(obj.direction))/(sigma^2));
                            gradient = voxelgrid.VEGFgrad(x,y,z);
                            if gradient == 0
                                gradient = [0 0 0];
                            end
                            if gradient ~= [0 0 0]
                                gradient = gradient./norm(gradient);
                            end
                            gradient_term = C_gradient.*(norm(gradient)/exp(1/sigma_g^2))*exp(dot(dir_vec./norm(dir_vec),gradient)/(sigma_g^2));
                         
                            if dot(dir_vec, obj.direction) < 0 %set points "behind" to neg infinity
                                coordinates{1,i,2}=0;
                                continue;
                            end
                            coordinates{1,i,2} = gamma*normcdf(tempVEGF - tipVEGF,0,sigma_VEGF)*persistence + gradient_term;
                            coordinates{1,i,1} = [x,y,z];
                            coordinates{1,i,3} = dir_vec;
   
                end
                
                score = [coordinates{1,:,2}];
                intervals = cumsum(score)./(max(cumsum(score)));
                maxindex = find(intervals >= rand, 1, 'first');
              
                %maxindices = find([coordinates{1,:,2}] == max([coordinates{1,:,2}]));
%                 if length(maxindices) == 1
                    dir_vec= coordinates{1,maxindex,3};

%                 else %choose randomly from maxes, if more than one. May need to add persistence bias!!
%                     %display('random!');
%                     random = randi(length(maxindices));
%                     new_coord = coordinates{1,maxindices(random),1};
%                     dir_vec= coordinates{1,maxindices(random),3};
% 
%                 end
            else
                %ADDED KERRI V2.2
                dir_vec = []; %There is no vector - search could not be completed 
                
            end
            
          end
          
          %% Canbranch method
          function [binary, binarycount] = canBranch(obj,XYZ, HypCell, pixelsize)
              global vegf_constant;
               if nargin < 3
                 error('Too few args')
               end
               %disp('!!!!!!!!!!CAN BRANCH!!!!!!!!!!!')
               %test if it is mature or not
               
                   % hypcell1 ={};
                   % hypcell2 = {};
                   %find cancer cell distances
                   segpos = obj.Node1./pixelsize;
                   segpos2 = obj.Node2./pixelsize;
                   %locations of hypoxic cells
                   XYZhyp  = XYZ(logical(HypCell),:)./pixelsize;
                   %Must have a hypoxic cell less than 250 microns away
                   % cdist1 = pdist2(segpos,20*XYZhyp);
                   % cdist2 = pdist2(segpos2,20*XYZhyp);
                   cdist1 = simpleDist(segpos,20*XYZhyp);
                   cdist2 = simpleDist(segpos2,20*XYZhyp);
                   % assert(isequal(cdist1(:),cdist1temp(:)))
                   % assert(isequal(cdist2(:),cdist2temp(:)))
                   hypcell1 = find(cdist1 <= 250,1);
                   hypcell2 = find(cdist2 <= 250,1);
                   if isempty(hypcell1) 
                       %don't branch 
                       VEGF1 = 0;
                   else
                       %branch
                       VEGF1 = vegf_constant;
                   end
                    if isempty(hypcell2) 
                       %don't branch 
                       VEGF2 = 0;
                   else
                       %branch
                       VEGF2 = vegf_constant;
                   end
             
                   
                   %depends on the cell's D114 status
                   switch obj.D114
                       case 0
                           probbranch =.2;%.06*(VEGF2.^3/(1+VEGF2.^3)); %(max(1-.5/VEGF2,0));
                           Vbranch = .5; % should be 1?
                       case 1
                           probbranch = .4;%.20;
                           Vbranch = 0;
                   end
                   binarycount = 0;
                   if (VEGF1 > Vbranch) || (VEGF2 > Vbranch) %Both nodes must have VEGF > .5 ng/ml
                       if obj.Branched == 0
                           binarycount = 1;
                       end
                       if (rand <= probbranch && obj.Branched == 0) && obj.Pepper ==0  %WILL ADDED 2/2
                           binary = 1;
%                            disp('branched')
%                            disp(obj.Node1)
%                            disp(obj.Node2)
%                            disp('capillary- ')
%                            disp(num2str(obj.CapillaryNum))
                       else
                           binary = 0;
                       end
                   else
                       binary = 0;
                   end
               

          end %canBranch
          
        
    end %End of Methods
    
end %End of segmentS

