function [newPs, skipIdx]=contourPropagate(Ps,sz, Options)

shrinkRate = Options.ShrinkPixelNum;
nPoints = Options.nPoints;
%lengthCanSkip = Options.lengthCanSkip;

    
skipIdx=[];

for i=1:1:numel(Ps)
    
%     if(Ps{i}.length<lengthCanSkip)
%         skipIdx=cat(1,skipIdx,i);
%         continue;
%     end
    
    P=Ps{i}.ctl;
    
    % Calculate distance between points
    dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];
    LastLength = dis(end);
    
    shrinkRate = Options.ShrinkPixelNum;
    if(LastLength<2*shrinkRate)
        if(isCloseToBoundary(P,sz(1),sz(2),Options.BoundThresh))
            skipIdx = cat(1,skipIdx,i);
            continue;
        else
            shrinkRate=1;
        end
    end
    
    if(isCloseToBoundary(P,sz(1),sz(2), Options.BoundThresh))
        LastLength = -1 * LastLength;
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
    
     % Calculate distance between points
    disK=cumsum(sqrt(sum((K(2:end,:)-K(1:end-1,:)).^2,2)));
    
    Ps{i}=struct('pts',K,'thickness',t,'length',disK(end),'targetLength',LastLength,...
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

