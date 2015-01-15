function [newPs, skipIdx]=contourPropagate(Ps,sz, Options)

shrinkRate = Options.ShrinkPixelNum;
nPoints = Options.nPoints;
lengthCanSkip = Options.lengthCanSkip;

    
skipIdx=[];

for i=1:1:numel(Ps)
    
    if(Ps{i}.length<lengthCanSkip)
        skipIdx=cat(1,skipIdx,i);
        continue;
    end
    
    P=Ps{i}.ctl;
    
    % Calculate distance between points
    dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];
    
    if(dis(end)<2*shrinkRate)
        skipIdx = cat(1,skipIdx,i);
        continue;
    end

    % Resample to make uniform points
    K=zeros(nPoints,2);
    K(:,1) = interp1(dis,P(:,1),linspace(shrinkRate,dis(end)-shrinkRate,nPoints));
    K(:,2) = interp1(dis,P(:,2),linspace(shrinkRate,dis(end)-shrinkRate,nPoints));
    %K(:,1) = interp1(dis,O(:,1),linspace(dis(end)*shrinkRate,dis(end)*(1-shrinkRate),nPoints));
    %K(:,2) = interp1(dis,O(:,2),linspace(dis(end)*shrinkRate,dis(end)*(1-shrinkRate),nPoints));
    
    % Clamp contour to boundary
    K(:,1)=min(max(K(:,1),1),sz(1));
    K(:,2)=min(max(K(:,2),1),sz(2));
    
    t = cellThickness(P,Ps{i}.seg,sz(1),sz(2));
    %len = Ps{i}.length;
    %len = Ps{i}.targetLength ;
    
    if(isCloseToBoundary(K,sz(1),sz(2), Options.BoundThresh))
        len = -1;
    else
        len = dis(end);
    end
    
    Ps{i}=struct('pts',K,'thickness',t,'length',0,'targetLength',len,...
    'strip1',[],'strip2',[],'region',[],'intensity',[],'normvec',[]);

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

