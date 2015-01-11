function [cellEachFrame,cMat]=updateCellEachFrame(cellEachFrame,newCellFrame,Ps,propagateIdx,tarMat,sz,divisionIDX)

cMat = zeros(sz);

%%% update confirmed cells %%%
sig=numel(newCellFrame);
for i=1:1:sig
    
    pp=newCellFrame{i}.ctl;
    idx=sub2ind(sz,pp(:,1),pp(:,2));
    cMat(idx)=i;
    
    pid=newCellFrame{i}.parent;
    
    if(isempty(pid)) % entering cell
        newCellFrame{i}=struct('length',newCellFrame{i}.length,'ctl',pp,'child',[],...
        'parent',[],'candi',[],'inflow',newCellFrame{i}.length,'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',newCellFrame{i}.seg);
    else
        newCellFrame{i}=struct('length',newCellFrame{i}.length,'ctl',pp,'child',[],...
        'parent',pid,'candi',[],'inflow',min([cellEachFrame{1}{pid}.length,...
        newCellFrame{i}.length]),'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',newCellFrame{i}.seg);
    end
   
end

clear idx pp i

%%% update evolved cells %%%
%skipIdx = [];
numProp = numel(propagateIdx);

for i=1:1:numProp
    % extract the centerline of region
    im=Ps{i}.region;
    [ctl, removedFlag] = pruneLine(bwmorph(im,'thin',Inf));
    
    if(removedFlag)
        %skipIdx = cat(2,skipIdx,i);
        continue;
    end
    
    ctlList=sortOneCellPixel(ctl);
    
    sig=sig+1;
    
    cMat(ctl>0)=sig;
    
    % associate the child of the corresponding cell in previous frame
    pid=propagateIdx(i);
    
    
    % insert the new cell
    mf = min([Ps{i}.length, cellEachFrame{1}{pid}.length]);
    
    tmp=struct('length',size(ctlList,1),'ctl',ctlList,'child',[],...
        'parent',pid,'candi',[],'inflow',mf,'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',im);
    newCellFrame=cat(2,newCellFrame,tmp);
    
    cellEachFrame{1}{pid}.child=sig;
    cellEachFrame{1}{pid}.outflow=mf;
end

if(~isempty(divisionIDX))
    if(numel(Ps)~=numProp + numel(divisionIDX))
        disp('error in divided cell contours')
        keyboard;
    end
    for i=numProp+1:1:numel(Ps)
        % extract the centerline of region
        im=Ps{i}.region;
        [ctl, removedFlag] = pruneLine(bwmorph(im,'thin',Inf));
        
        if(removedFlag)
            %skipIdx = cat(2,skipIdx,i);
            continue;
        end
        
        ctlList=sortOneCellPixel(ctl);
        
        sig=sig+1;
        
        cMat(ctl>0)=sig;
        
        % associate the child of the corresponding cell in previous frame
        pid=propagateIdx(divisionIDX(i-numProp));

        % insert the new cell
        mf = Ps{i}.length;
        
        tmp=struct('length',size(ctlList,1),'ctl',ctlList,'child',[],...
            'parent',pid,'candi',[],'inflow',mf,'outflow',0,'relaxinCost',0,...
            'relaxOutCost',0,'seg',im);
        newCellFrame=cat(2,newCellFrame,tmp);
        
        cellEachFrame{1}{pid}.child=cat(2,sig,cellEachFrame{1}{pid}.child);
        cellEachFrame{1}{pid}.outflow=min([cellEachFrame{1}{pid}.length, mf+cellEachFrame{1}{pid}.outflow]);
    end
end

%%% update matching between 2 and 3
for i=1:1:numel(cellEachFrame{3})
    cellEachFrame{3}{i}.inflow=0;
    cellEachFrame{3}{i}.parent=[];
end


[srcCellList,tarCellList]=local_EMD(newCellFrame,cellEachFrame{3}, cMat, tarMat);
cellEachFrame{2}=srcCellList;
cellEachFrame{3}=tarCellList;
