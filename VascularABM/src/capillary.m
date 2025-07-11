%Class capillary, is an easy way to keep track of separate capillaries
%Version 2.1
%Last user: Will
%Date: 1.24

classdef capillary<handle
   properties
       SegmentList = cell(100,1); %List of segments
       ActiveList = zeros(100,1); %List of active segments
       numSeg = 0; %Number of segments
       BranchList = zeros(100,1); %List of branched segments
       %Listpos = 0;
   end
   
   methods
       function [obj] = addSegment(obj,newSegment)
          emptylist = cellfun(@(A) any(isempty(A(:))), obj.SegmentList);
          emptypos = find (emptylist, 1, 'first');
          %newSegment.CapillaryNum = obj.Listpos;
          obj.SegmentList{emptypos} = newSegment;
          obj.numSeg = obj.numSeg + 1;
       end
       
       function [obj] = deactivateCapillary(obj)
          for i = 1:length(obj.SegmentList(~cellfun('isempty',obj.SegmentList)))
                obj.SegmentList{i}.Pepper = 1;
          end
       end
       
   end
    
end