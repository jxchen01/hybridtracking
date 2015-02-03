function [newPs, skipIdx]=contourPropagate(cellEachFrame, idxList,I, Options)

sz=size(I);
nPoints = Options.nPoints;
Ps=cellEachFrame{1}(idxList);

skipIdx=[];

for i=1:1:numel(Ps)
    
    P=Ps{i}.ctl; % pixel-level accuracy (all connected grid points)
    
    % keep record of last cell and decide shrinkrate
    dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];
    LastLength = dis(end);
    t = cellThickness(P,Ps{i}.seg,sz(1),sz(2));
    P0=P;
    
    shrinkRate = Options.ShrinkPixelNum;
    if(LastLength<=2*shrinkRate+3)
        if(isCloseToBoundary(P,sz(1),sz(2),Options.BoundThresh))
            skipIdx = cat(1,skipIdx,i);
            continue;
        else
            shrinkRate=max([1,floor(LastLength*0.35)]);
        end
    end
    
    if(isCloseToBoundary(P,sz(1),sz(2), Options.BoundThresh))
        LastLength = -1 * LastLength;
    elseif(isfield(Ps{i},'dangerLength'))
        LastLength = Ps{i}.dangerLength;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% calculate the moving direction %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cid = Ps{i}.child;
    if(numel(cid)>0)
        % get corresponding cells in future frames
        cellNext=[];
        cellFuture=[];
        cf=Ps{i}.cumFlow;
        for kk=1:1:numel(cid)
            if(abs(cid(kk) - round(cid(kk)))<1e-8)
                cellNext=cat(2,cellNext,cid(kk));
            else
                cellFuture = cat(1,cellFuture,[floor(cid(kk)), uint16((cid(kk) - floor(cid(kk)))*1000), cf(kk)]);
            end
        end
        
        if(~isempty(cellNext))
            if(numel(cellNext)==1)
                P=cellEachFrame{2}{cellNext}.ctl;
            else
                P=mergeCells(cellEachFrame{2}(cellNext));
            end
        else
            [mf,midx]=max(cellFuture(:,3));
            P = cellMorphing(P, cellEachFrame{cellFuture(midx,1)}{cellFuture(midx,2)}.ctl,...
                mf, cellFuture(midx,1)-1);
            clear mf midx
        end
        
        clear cellNext cellFuture cf kk
    end
    
    %%%%%%
    % Note: P has been updated
    %%%%%%
    
    % Calculate distance between points
    dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];

    % Resample to make uniform points
    K=zeros(nPoints,2);
    K(:,1) = interp1(dis,P(:,1),linspace(shrinkRate,dis(end)-shrinkRate,nPoints));
    K(:,2) = interp1(dis,P(:,2),linspace(shrinkRate,dis(end)-shrinkRate,nPoints));
    
    % Clamp contour to boundary
    K(:,1)=min(max(K(:,1),1),sz(1));
    K(:,2)=min(max(K(:,2),1),sz(2));
    
     % Calculate distance between points
    disK=cumsum(sqrt(sum((K(2:end,:)-K(1:end-1,:)).^2,2)));
    
    % copy current information so that the evolved contour can be examined
    SingleContour=false(sz);
    [R1,R2,NV]=getRibbon(K,t,sz(1),sz(2));
    contourList=[R1(:,:);R2(end:-1:1,:);R1(1,:)];
    for ic=2:1:size(contourList,1)
        [xp,yp]=bresenham(contourList(ic,1),contourList(ic,2),contourList(ic-1,1),contourList(ic-1,2));
        pidx=sub2ind(sz,xp,yp);
        SingleContour(pidx)=true;
    end
    SingleContour=imfill(SingleContour,'holes');
    interiorIntensity = mean(I(SingleContour>0));

    Ps{i}=struct('pts',K,'thickness',t,'length',disK(end),'targetLength',LastLength,...
    'strip1',R1,'strip2',R2,'region',SingleContour,'intensity',interiorIntensity,...
    'normvec',NV,'LastFrameIntensity',interiorIntensity,'LastFramePts',P0);
    % Note: 'LastFrameIntensity' is actually the intensity before evolution
    % in current frame

    clear O len K P dis t
end

if(~isempty(skipIdx))
    numCell=numel(Ps)-numel(skipIdx);
    newPs=cell(1,numCell);
    idx = setdiff(1:1:numel(Ps), skipIdx);
    for i=1:1:numCell
        newPs{i}=Ps{idx(i)};
    end
else
    newPs = Ps;
end

