function [cellEachFrame,cMat,maxID]=updateCellEachFrame(cellEachFrame,oldCellFrame,Ps,propagateIdx,tarMat,sz,divisionIDX,maxID,Options)

cMat = zeros(sz);
newCellFrame=cell(1,0);
%%% update evolved cells %%%
numProp = numel(propagateIdx);
evolvedMap = zeros(sz);
sig=0;
for i=1:1:numProp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check whether the contour should be removed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Case 1: If it is shorter than a threshold, close to boundary,
    %         also becoming shorter than its length in the last frame
    if((Ps{i}.length < Options.lengthCanSkip)  ...
            && isCloseToBoundary(Ps{i}.pts,sz(1),sz(2),Options.BoundThresh)...
            && Ps{i}.length < 0.5+abs(Ps{i}.targetLength))
        continue;
    end
    
    % Case 2: If the interior intensity differ too much from that 
    %         in the last frame. (Evolution Error)
    if(CellDist(Ps{i})>2.5 || Ps{i}.intensity/(Ps{i}.LastFrameIntensity+0.000001)< 0.7)
        disp('check evolution');
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract the segmentation region
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im=Ps{i}.region;
    
    %%% add two heads %%%
    tmpHead=zeros(sz);
    pts = Ps{i}.pts;
    x1=round(pts(1,1));y1=round(pts(1,2));
    if(x1<1), x1=1; elseif(x1>sz(1)), x1=sz(1); end
    if(y1<1), y1=1; elseif(y1>sz(2)), y1=sz(2); end
    tmpHead(x1,y1)=1;
    x2=round(pts(end,1));y2=round(pts(end,2));
    if(x2<1), x2=1; elseif(x2>sz(1)), x2=sz(1); end
    if(y2<1), y2=1; elseif(y2>sz(2)), y2=sz(2); end
    tmpHead(x2,y2)=1;
    
    se= strel('disk',double(max([1,round(Ps{i}.thickness)])),0);
    im = im | imdilate(tmpHead,se);
    
    evolvedMap = evolvedMap | im;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert contour control points to centerline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pts= Ps{i}.pts;
    ctl=zeros(sz);
    for pid=1:1:size(pts,1)-1
        [px,py]=bresenham(pts(pid,1),pts(pid,2),pts(pid+1,1),pts(pid+1,2));
        pidx = sub2ind(sz,px,py);
        ctl(pidx)=1;
    end
    ctl = bwmorph(ctl,'thin',Inf);
    ctlList=sortOneCellPixel(ctl);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add the centerline the cMat
    sig=sig+1;
    cMat(ctl>0)=sig;
    
    % associate the child of the corresponding cell in previous frame
    pid=propagateIdx(i);  
    
    % compute the flow amount
    mf = min([Ps{i}.length, cellEachFrame{1}{pid}.length]);
    
    % build the new cell
    tmp=struct('length',size(ctlList,1),'ctl',ctlList,'child',[],...
        'parent',pid,'candi',[],'inflow',mf,'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',im, 'id',cellEachFrame{1}{pid}.id);
    
    % when the length decreases too much, fire alwarm by setting dangerLength
    if(Ps{i}.length<max([0.8*Ps{i}.targetLength, Ps{i}.targetLength-4])...
            && ~isCloseToBoundary(Ps{i}.pts,sz(1),sz(2),Options.BoundThresh))
        tmp.dangerLength=Ps{i}.targetLength;
    end
        
    % insert the new cell
    newCellFrame=cat(2,newCellFrame,tmp);
    
    % update the parents the new cell in previous frame
    cellEachFrame{1}{pid}.child=sig;
    cellEachFrame{1}{pid}.outflow=mf;
end

if(~isempty(divisionIDX))
    if(numel(Ps)~=numProp + numel(divisionIDX))
        disp('error in divided cell contours')
        keyboard;
    end
    for i=numProp+1:1:numel(Ps)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check whether the contour should be removed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Case 1: If it is shorter than a threshold, close to boundary,
        %         also becoming shorter than its length in the last frame
        if((Ps{i}.length < Options.lengthCanSkip)  ...
                && isCloseToBoundary(Ps{i}.pts,sz(1),sz(2),Options.BoundThresh)...
                && Ps{i}.length < 0.5+abs(Ps{i}.targetLength))
            continue;
        end
        
        % Case 2: If the interior intensity differ too much from that
        %         in the last frame. (Evolution Error)
        if(CellDist(Ps{i})>2.5 || Ps{i}.intensity/(Ps{i}.LastFrameIntensity+0.000001)< 0.7)
            disp('check evolution');
            continue;
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract the segmentation region
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        im=Ps{i}.region;
        
        %%% add two heads %%%
        tmpHead=zeros(sz);
        pts = Ps{i}.pts;
        x1=round(pts(1,1));y1=round(pts(1,2));
        if(x1<1), x1=1; elseif(x1>sz(1)), x1=sz(1); end
        if(y1<1), y1=1; elseif(y1>sz(2)), y1=sz(2); end
        tmpHead(x1,y1)=1;
        x2=round(pts(end,1));y2=round(pts(end,2));
        if(x2<1), x2=1; elseif(x2>sz(1)), x2=sz(1); end
        if(y2<1), y2=1; elseif(y2>sz(2)), y2=sz(2); end
        tmpHead(x2,y2)=1;
        
        se= strel('disk',double(max([1,round(Ps{i}.thickness)])),0);
        im = im | imdilate(tmpHead,se);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert contour control points to centerline
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pts= Ps{i}.pts;
        ctl=zeros(sz);
        for pid=1:1:size(pts,1)-1
            [px,py]=bresenham(pts(pid,1),pts(pid,2),pts(pid+1,1),pts(pid+1,2));
            pidx = sub2ind(sz,px,py);
            ctl(pidx)=1;
        end
        ctl = bwmorph(ctl,'thin',Inf);
        ctlList=sortOneCellPixel(ctl);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add the centerline the cMat
        sig=sig+1;    
        cMat(ctl>0)=sig;
        
        % associate the child of the corresponding cell in previous frame
        pid=propagateIdx(divisionIDX(i-numProp));

        % compute the flow amount
        mf = Ps{i}.length;
        
        % build the new cell
        tmp=struct('length',size(ctlList,1),'ctl',ctlList,'child',[],...
            'parent',pid,'candi',[],'inflow',mf,'outflow',0,'relaxinCost',0,...
            'relaxOutCost',0,'seg',im, 'id',cellEachFrame{1}{pid}.id);
        
        % insert the new cell
        newCellFrame=cat(2,newCellFrame,tmp);
        
        cellEachFrame{1}{pid}.child=cat(2,sig,cellEachFrame{1}{pid}.child);
        cellEachFrame{1}{pid}.outflow=min([cellEachFrame{1}{pid}.length, mf+cellEachFrame{1}{pid}.outflow]);
    end
end

%%% update confirmed cells %%%
for i=1:1:numel(oldCellFrame)
    if(~isfield(oldCellFrame{i},'case'))
        keyboard;
    end
    if(oldCellFrame{i}.case<0) % re-appearing cell
        tmpSeg= oldCellFrame{i}.seg;
        if(nnz(tmpSeg & evolvedMap)>3)
            continue;
        end
    end
    
    sig=sig+1;
    
    pp=oldCellFrame{i}.ctl;
    idx=sub2ind(sz,pp(:,1),pp(:,2));
    cMat(idx)=sig;
    
    pid=oldCellFrame{i}.parent;
    
    if(isempty(pid)) % entering cell or re-appearing
        maxID = maxID + 1;
        tmp=struct('length',oldCellFrame{i}.length,'ctl',pp,'child',[],...
        'parent',[],'candi',[],'inflow',oldCellFrame{i}.length,'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',oldCellFrame{i}.seg,'id',maxID);
    else
        tmp=struct('length',oldCellFrame{i}.length,'ctl',pp,'child',[],...
        'parent',pid,'candi',[],'inflow',min([cellEachFrame{1}{pid}.length,...
        oldCellFrame{i}.length]),'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',oldCellFrame{i}.seg, 'id', cellEachFrame{1}{pid}.id);
    end
    
    newCellFrame = cat(2, newCellFrame,tmp);   
end

clear idx pp i


%%% update matching between 2 and 3
for i=1:1:numel(cellEachFrame{3})
    cellEachFrame{3}{i}.inflow=0;
    cellEachFrame{3}{i}.parent=[];
end

[srcCellList,tarCellList]=local_EMD(newCellFrame,cellEachFrame{3}, cMat, tarMat, Options);
cellEachFrame{2}=srcCellList;
cellEachFrame{3}=tarCellList;
