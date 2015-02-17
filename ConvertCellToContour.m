function [Ps, newCellFrame, BMap, idxList]=ConvertCellToContour(cellEachFrame,I, Options)

sz=size(I);
cFrame = cellEachFrame{2};
pFrame = cellEachFrame{1};

newCellFrame = cell(1,0);
propagateIdx =[];

BMap = zeros(sz);

for i=1:1:numel(cFrame);
    
    if(isempty(cFrame{i}.parent)) 
        
        %%%% a short cell will not be considered for entering %%%%
        %%%% short cell may introduce big trouble in evolution %%%%%
        %%%% because control points may lump together %%%%%
        if(cFrame{i}.length < Options.pixelCanSkip)
            continue;
        end
       
        if(isCloseToBoundary(cFrame{i}.ctl,sz(1),sz(2), Options.BoundThresh))
            if( confirmEntry(cellEachFrame(1,2:1:end),i,Options.lengthCanSkip,1))  
                % confirmed entering cell
                newCellFrame = cat(2,newCellFrame, cFrame{i});
                newCellFrame{end}.case = 2; % entering
                
%                 % only confirmed entering cell can contribute to the BMap,
%                 % besides good segmentation
%                 BMap = BMap | cFrame{i}.seg;
            end
        elseif(confirmEntry(cellEachFrame(1,2:1:end),i,Options.lengthCanSkip,2))          
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
end

idxList = setdiff(1:1:numel(pFrame), propagateIdx);
[Ps, skipIdx] = contourPropagate(cellEachFrame,idxList,I, Options);
if(~isempty(skipIdx))
    idxList = setdiff( idxList, idxList(skipIdx));
end


