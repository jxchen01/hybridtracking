function [Ps, newCellFrame, BMap, idxList]=ConvertCellToContour(cellEachFrame,sz)

cFrame = cellEachFrame{2};
pFrame = cellEachFrame{1};

newCellFrame = cell(1,0);
propagateIdx =[];

BMap = zeros(sz);

for i=1:1:numel(cFrame);
    
    if(isempty(cFrame{i}.parent)) 
        if(isCloseToBoundary(cFrame{i}.ctl,sz(1),sz(2)) && confirmEntry(cellEachFrame(1,2:1:end),i))
            % confirmed entering cell
            newCellFrame = cat(2,newCellFrame, cFrame{i});
            BMap = BMap | cFrame{i}.seg;
        end
    else
        pidx = cFrame{i}.parent;
        if(numel(pidx)==1  && abs(pFrame{pidx}.length - ...
                cFrame{i}.length)<= max(4,0.08*pFrame{pidx}.length))
            % confirmed good segmentation
            newCellFrame = cat(2, newCellFrame, cFrame{i});
            propagateIdx = cat(2,propagateIdx, pidx);
            
            BMap = BMap | cFrame{i}.seg;
        end       
    end
    
end

idxList = setdiff(1:1:numel(pFrame), propagateIdx);
[Ps, skipIdx] = contourPropagate(pFrame(idxList),10,sz);
if(~isempty(skipIdx))
    idxList = setdiff( idxList, idxList(skipIdx));
end


