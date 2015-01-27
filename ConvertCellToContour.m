function [Ps, newCellFrame, BMap, idxList]=ConvertCellToContour(cellEachFrame,I, Options)

sz=size(I);
cFrame = cellEachFrame{2};
pFrame = cellEachFrame{1};

newCellFrame = cell(1,0);
propagateIdx =[];

BMap = zeros(sz);

for i=1:1:numel(cFrame);
    
    if(isempty(cFrame{i}.parent)) 
       
        if(isCloseToBoundary(cFrame{i}.ctl,sz(1),sz(2), Options.BoundThresh)...
                && confirmEntry(cellEachFrame(1,2:1:end),i,1))
            
            % confirmed entering cell 
            newCellFrame = cat(2,newCellFrame, cFrame{i});
            newCellFrame{end}.case = 2; % entering
            
            % only confirmed entering cell can contribute to the BMap
            BMap = BMap | cFrame{i}.seg; 
        elseif(confirmEntry(cellEachFrame(1,2:1:end),i,2))          
            % confirmed re-appearing cell 
            newCellFrame = cat(2,newCellFrame, cFrame{i});
            newCellFrame{end}.case = -1; % re-appearing         
        end
    else
        pidx = cFrame{i}.parent;
        if(numel(pidx)==1  && abs(pFrame{pidx}.length - ...
                cFrame{i}.length)<= max(4,0.08*pFrame{pidx}.length))
            % confirmed good segmentation
            newCellFrame = cat(2, newCellFrame, cFrame{i});
            newCellFrame{end}.case =1 ; % good
            propagateIdx = cat(2,propagateIdx, pidx);
            
            BMap = BMap | cFrame{i}.seg;
        end       
    end
%     if(numel(newCellFrame)==82)
%         keyboard;
%     end
end

idxList = setdiff(1:1:numel(pFrame), propagateIdx);
[Ps, skipIdx] = contourPropagate(pFrame(idxList),I, Options);
if(~isempty(skipIdx))
    idxList = setdiff( idxList, idxList(skipIdx));
end


