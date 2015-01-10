function [cellEachFrame,cMat]=updateCellEachFrame(cellEachFrame,newCellFrame,Ps,propagateIdx,tarMat,sz)

cMat = zeros(sz);

%%% update confirmed cells %%%
sig=numel(newCellFrame);
for i=1:1:sig
    
    pp=newCellFrame{i}.ctl;
    idx=sub2ind(sz,pp(:,1),pp(:,2));
    cMat(idx)=i;
    
    pid=newCellFrame{i}.parent;
    
    newCellFrame{i}=struct('length',newCellFrame{i}.length,'ctl',pp,'child',[],...
        'parent',pid,'candi',[],'inflow',min([cellEachFrame{1}{pid}.length,...
        newCellFrame{i}.length]),'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',newCellFrame{i}.seg);
end

clear idx pp i

%%% update evolved cells %%%
for i=1:1:numel(Ps)
    sig=sig+1;
    
    pid=propagateIdx(i);
    cellEachFrame{1}{pid}.child=sig;
    
    im=Ps{i}.region;
    ctl = pruneLine(bwmorph(im,'thin',Inf));
    ctlList=sortOneCellPixel(ctl);
    
    mf = min([Ps{i}.length, cellEachFrame{1}{pid}.length]);
    
    tmp=struct('length',size(ctlList,1),'ctl',ctlList,'child',[],...
        'parent',pid,'candi',[],'inflow',mf,'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',im);
    newCellFrame=cat(2,newCellFrame,tmp);
    
    cellEachFrame{1}{pid}.outflow=mf;
    
    cMat(ctl>0)=sig;
end

%%% update matching between 2 and 3
for i=1:1:numel(cellEachFrame{3})
    cellEachFrame{3}{i}.inflow=0;
    cellEachFrame{3}{i}.parent=[];
end


[srcCellList,tarCellList]=local_EMD(newCellFrame,cellEachFrame{3}, cMat, tarMat);
cellEachFrame{2}=srcCellList;
cellEachFrame{3}=tarCellList;
