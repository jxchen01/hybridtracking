function [newPs, skipIdx]=contourPropagate(cellEachFrame, idxList,I, Options)

sz=size(I);
nPoints = Options.nPoints;
Ps=cellEachFrame{1}(idxList);
Ps0=Ps;

skipIdx=[];

for i=1:1:numel(Ps)

    P=Ps{i}.ctl; % pixel-level accuracy (all connected grid points)
    
    % propagate target length 
    dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];
    LastLength = dis(end);
    if(isfield(Ps{i},'dangerLength') && Ps{i}.dangerLength>1e-5)
        LastLength = Ps{i}.dangerLength;
    end
    
    if(isCloseToBoundary(P,sz(1),sz(2), Options.BoundThresh))
        LastLength = -1 * LastLength;
    end
    
    % propagate thickness
    t = cellThickness(P,Ps{i}.seg,sz(1),sz(2));
    
    P0=P;   
    
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
                cellNext=cat(1,cellNext,[cid(kk),cf(kk)]);
            else
                cellFuture = cat(1,cellFuture,[floor(cid(kk)), uint16((cid(kk) - floor(cid(kk)))*1000), cf(kk)]);
            end
        end
        
        % apply different methods to do morphing in different cases
        if(~isempty(cellNext)) 
            % do morphing directly with the cells in the next frame
            if(size(cellNext,1)==1)  
                if(numel(cellEachFrame{2}{cellNext(1,1)}.parent)==1)
                    % 1-to-1 matching
                    if(cellNext(1,2)<Options.bodyRatio*cellEachFrame{2}{cellNext(1,1)}.length)
                        % match to a larger one
                        P=cellMorphing(P, cellEachFrame{2}{cellNext(1,1)}.ctl,cellNext(1,2), 1,sz);
                    else
                        % match to a smaller one or of similar size
                        P=cellEachFrame{2}{cellNext(1,1)}.ctl;
                    end  
                else
                    % N-to-1 matching
                    pList=cellEachFrame{2}{cellNext(1,1)}.parent;
                    pnum=numel(pList);
                    pidx=[];
                    pflag=0;
                    for kk=1:1:pnum
                        a=find(idxList==pList(kk));
                        if(~isempty(a))
                            pidx=cat(2,pidx,a);
                            if(a==i)
                                pflag=numel(pidx);
                            end
                        end
                    end
                    if(numel(pidx)>1)
                        P=multiMorphing(Ps0(pidx),cellEachFrame{2}{cellNext(1,1)},pflag,sz);
                        if(isempty(P))
                            skipIdx = cat(1,skipIdx,i);
                            continue;
                        end
                    else % some of the parents may have been confirmed
                         % then, they should not be included for morphing
                        P=cellMorphing(P, cellEachFrame{2}{cellNext(1,1)}.ctl,cellNext(1,2), 1,sz);
                    end
                end
            else
                % 1-to-N matching
                P=mergeCells(cellEachFrame{2}(cellNext(:,1)));
            end
        else
            [mf,midx]=max(cellFuture(:,3));
            P = cellMorphing(P, cellEachFrame{cellFuture(midx,1)}{cellFuture(midx,2)}.ctl,...
                mf, cellFuture(midx,1)-1, sz);
            clear mf midx
        end
        
        clear cellNext cellFuture cf kk
        
        %%%%% remove, if the initial position has too few points %%%%%
        %%%% this could be caused by small amount of flow between cells %%%
        if(size(P,1)<5)
            skipIdx = cat(1,skipIdx,i);
            continue;
        end
        
        % Calculate distance between points
        dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];
        
        % determine shrinkRate
        shrinkRate = Options.ShrinkPixelNum;
        shrinkLengthRation = 0.3;
    else   
        %%%% remove leaving cells %%%%%%
        if(isCloseToBoundary(P,sz(1),sz(2),Options.BoundThresh))
            %if(LastLength<=Options.leavingLength && isCloseToBoundary(P,sz(1),sz(2),Options.BoundThresh))
            skipIdx = cat(1,skipIdx,i);
            continue;
        end
        
        % determine shrinkRate
        shrinkRate = Options.ShrinkPixelNum + 12;
        shrinkLengthRation = 0.4;
    end
    
    %%%%%%
    % Note: P and dis has been updated
    %%%%%%
    
    %%%%% remove short cells %%%%%%
    if(dis(end)<Options.lengthCanSkip)
        skipIdx = cat(1,skipIdx,i);
        continue;
    end
    
    if(dis(end)<2*shrinkRate+5)
        shrinkRate=max([1,floor(dis(end)*shrinkLengthRation)]);
    end
    
    % Resample to make uniform points
    K=zeros(nPoints,2);
    
    %%%%%%%%%%%%%%
    % if close to boundary, shrink one side; otherwise, shrink both sides %
    %%%%%%%%%%%%%%
    if(isCloseToBoundary(P,sz(1),sz(2), Options.BoundThresh))
        % determine which side is closer to boundary
        if(distToBoundary(P(1,1),P(1,2),sz(1),sz(2)) < distToBoundary(P(end,1),P(end,2),sz(1),sz(2)))
            K(:,1) = interp1(dis,P(:,1),linspace(0,dis(end)-shrinkRate,nPoints));
            K(:,2) = interp1(dis,P(:,2),linspace(0,dis(end)-shrinkRate,nPoints));
        else
            K(:,1) = interp1(dis,P(:,1),linspace(shrinkRate,dis(end),nPoints));
            K(:,2) = interp1(dis,P(:,2),linspace(shrinkRate,dis(end),nPoints));
        end
    else
        K(:,1) = interp1(dis,P(:,1),linspace(shrinkRate,dis(end)-shrinkRate,nPoints));
        K(:,2) = interp1(dis,P(:,2),linspace(shrinkRate,dis(end)-shrinkRate,nPoints));       
    end  
    
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

end



function dist = distToBoundary(px,py,xdim,ydim)
    hx=min([px,xdim-px]);
    hy=min([py,ydim-py]);
    
    dist = min([hx,hy]);
end

