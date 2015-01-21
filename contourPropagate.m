function [newPs, skipIdx]=contourPropagate(Ps,I, Options)

sz=size(I);
nPoints = Options.nPoints;
%lengthCanSkip = Options.lengthCanSkip;

skipIdx=[];

for i=1:1:numel(Ps)
    
    P=Ps{i}.ctl;
    
    % Calculate distance between points
    dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];
    LastLength = dis(end);
    
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
    
    % compute cell thichness
    t = cellThickness(P,Ps{i}.seg,sz(1),sz(2));

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
    'normvec',NV,'LastFrameIntensity',interiorIntensity,'LastFramePts',P);

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

