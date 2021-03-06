function cellFrame = Global_EMD(cellFrame, matFrame, algOptions)

%%%%%%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numFrame = numel(matFrame);
[dimx, dimy] = size(matFrame{1}.Mat);

bodyRatio = algOptions.bodyRatio;
searchDepth = algOptions.searchDepth;
candiRadiusK = algOptions.candiRadiusK;

se=cell(1,searchDepth-1);
for i=1:1:searchDepth-1
    se{i}=strel('disk',candiRadiusK(i));
end

%options = optimset('Algorithm','simplex','Display', 'off', 'Diagnostics',...
%    'off','LargeScale', 'off', 'Simplex', 'on');
options = optimset('Display', 'off');

tarMat=zeros(dimx,dimy,searchDepth-1);
flowInfo = cell(500);
numFlowEdge = 0;
flowConst = [];
outNum = 0;
flowCapIdx = zeros(numFrame, 1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% check outflow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sout=0;
%%%%%%%%%%%%%%%%%%%%%% loop through all frames %%%%%%%%%%%%%%%%%%%%%%%
for frameK = 1:1:numFrame-1
    srcCellFrame = cellFrame{frameK};
    tmpCellNum = numel(srcCellFrame);
    srcMat=matFrame{frameK}.Mat;
    for i=1:1:searchDepth-1
        if(i+frameK+1>numFrame)
            break;
        end
        tmp=matFrame{i+frameK+1}.Mat;
        tarMat(:,:,i)=tmp(:,:);
    end
    clear tmp 
    
    for i=1:1:tmpCellNum
        tmpCell = srcCellFrame{i};
        if(tmpCell.outflow >= bodyRatio*tmpCell.length)
            continue;
        end
        
        % If reaching here, it means the node has not been satisfied.
        flowConst = cat(1,flowConst,(tmpCell.length-tmpCell.outflow));
        outNum = outNum+1;
        Sout=Sout+(tmpCell.length-tmpCell.outflow);
          
        % search in layer K+1
        tmpCandi=setdiff(tmpCell.candi,tmpCell.child);
        for j=1:1:numel(tmpCandi)
            if(cellFrame{frameK+1}{tmpCandi(j)}.inflow >= bodyRatio*cellFrame{frameK+1}{tmpCandi(j)}.length)
                continue;
            end
            
            flowCap = min([tmpCell.length-tmpCell.outflow , ...
                cellFrame{frameK+1}{tmpCandi(j)}.length-cellFrame{frameK+1}{tmpCandi(j)}.inflow]);
            
            % if reach here, it means a potential flow is found
            tmpD = ComputeDist(tmpCell.ctl,tmpCell.length, cellFrame{frameK+1}{tmpCandi(j)}.ctl,cellFrame{frameK+1}{tmpCandi(j)}.length,0);
            if(tmpD<500)
                numFlowEdge = numFlowEdge + 1;
                flowInfo{numFlowEdge} = struct('head',[frameK,i],...
                    'tail',[frameK+1,tmpCandi(j)],'edgeCost',tmpD,...
                    'ub',flowCap,'outConstIdx', outNum);
            end
        end
        
        % search in layer k+2 ... 
        for kk=2:1:searchDepth
            idxFM=frameK+kk;
            if(idxFM>numFrame)
                break;
            end
            tar_idx=kk-1;
            tmpSingle=(srcMat==i);
            tmpSingle=imdilate(tmpSingle,se{tar_idx});
            tmpAND = tmpSingle.*tarMat(:,:,tar_idx);
            tmpCandi = unique(nonzeros(tmpAND));
            
            for j=1:1:length(tmpCandi)
                if(cellFrame{idxFM}{tmpCandi(j)}.inflow >= bodyRatio*cellFrame{idxFM}{tmpCandi(j)}.length )
                    continue;
                end
                
                 % if reach here, it means a potential flow is found
                flowCap = min([tmpCell.length-tmpCell.outflow , cellFrame{idxFM}{tmpCandi(j)}.length-cellFrame{idxFM}{tmpCandi(j)}.inflow]);                

                tmpD = ComputeDist(tmpCell.ctl,tmpCell.length, cellFrame{idxFM}{tmpCandi(j)}.ctl,cellFrame{idxFM}{tmpCandi(j)}.length,0);
                if(tmpD<500)
                    numFlowEdge = numFlowEdge + 1;
                    flowInfo{numFlowEdge} = struct('head',[frameK,i], ...
                        'tail',[idxFM,tmpCandi(j)],'edgeCost',tmpD,...
                        'ub',flowCap,'outConstIdx', outNum);
                end
            end
        end
        
        % add a flow to sink node
        tmpCost=tmpCell.relaxOutCost;
        if(tmpCost<0)
            tmpCost=-tmpCost;
        elseif(tmpCost==0)
            keyboard
        end
        numFlowEdge = numFlowEdge + 1;
        flowInfo{numFlowEdge} = struct('head',[frameK,i], 'tail', [numFrame+1, 0],...
            'edgeCost',tmpCost,'ub',tmpCell.length-tmpCell.outflow,'outConstIdx', outNum);      
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% check inflow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inNum = outNum;
Sinf=0;

for frameK=2:numFrame
    srcCellFrame = cellFrame{frameK};
    for i=1:1:numel(srcCellFrame)
        tmpCell = srcCellFrame{i};
        if(tmpCell.inflow >= bodyRatio*tmpCell.length )
            continue;
        end
        
        % add a flow from source node
        tmpCost=tmpCell.relaxInCost;
        if(tmpCost<0)
            tmpCost=-tmpCost;
        end
        
        numFlowEdge = numFlowEdge + 1;
        flowInfo{numFlowEdge} = struct('head',[0,0], 'tail', [frameK, i],...
            'edgeCost',tmpCost,'ub',tmpCell.length-tmpCell.inflow,'outConstIdx',0);
        
        flowConst = cat(1,flowConst,(tmpCell.length-tmpCell.inflow));
        Sinf=Sinf+(tmpCell.length-tmpCell.inflow);
        inNum = inNum+1;
        flowCapIdx(frameK,i)=inNum;
    end
    
end

clear tmpCost
%%%%%%%%%%%%%%%%%%%%%%% build the LP problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight=zeros(1,numFlowEdge);
lb=zeros(1,numFlowEdge);
ub=zeros(1,numFlowEdge);
Aneq=zeros(inNum,numFlowEdge);
Aeq=ones(1,numFlowEdge);
beq=min([Sinf,Sout]);

for i=1:1:numFlowEdge
    tmpEdge=flowInfo{i};
    % update weight
    weight(i)=tmpEdge.edgeCost;
    
    % update Aneq
    if(tmpEdge.outConstIdx>0)
        Aneq(tmpEdge.outConstIdx,i)=1;
        
        tmpTail=tmpEdge.tail;
        if(tmpTail(1)<1+numFrame)
            tmpIdx=flowCapIdx(tmpTail(1),tmpTail(2));
            Aneq(tmpIdx,i)=1;     
        end
    else
        tmpTail=tmpEdge.tail;
        tmpIdx=flowCapIdx(tmpTail(1),tmpTail(2));
        Aneq(tmpIdx,i)=1;
    end
    
    % update ub
    ub(i)=tmpEdge.ub;
end
clear tmpEdge tmpIdx tmpTail 

[xval,~,exitflag,output] = linprog(weight,Aneq,flowConst,Aeq,beq,lb,[],[],options);

if(exitflag~=1)
    disp(output.message);
    error('error in optimization');
end

for i=1:1:numFlowEdge
    if(xval(i)>0.1)       
        tmpHead=flowInfo{i}.head;
        tmpTail=flowInfo{i}.tail;          
        if(tmpHead(1)~=0 && tmpTail(1)~=numFrame+1)   
            check1= xval(i)>bodyRatio*(cellFrame{tmpHead(1)}{tmpHead(2)}.length...
                -cellFrame{tmpHead(1)}{tmpHead(2)}.outflow);
            check2= xval(i)>bodyRatio*(cellFrame{tmpTail(1)}{tmpTail(2)}.length-...
                cellFrame{tmpTail(1)}{tmpTail(2)}.inflow);
            check3= xval(i)>0.25*cellFrame{tmpHead(1)}{tmpHead(2)}.length;
            check4= xval(i)>0.25*cellFrame{tmpTail(1)}{tmpTail(2)}.length;
            if( (check1||check2) && (check3||check4) )
       
                tmp=cellFrame{tmpHead(1)}{tmpHead(2)}.child;
                tmp2=cellFrame{tmpTail(1)}{tmpTail(2)}.parent;
                if(tmpTail(1)==tmpHead(1)+1)
                    tmp=cat(2,tmp,tmpTail(2));
                    tmp2=cat(2,tmp2,tmpHead(2));
                else
                    tmp=cat(2,tmp,tmpTail(1)+tmpTail(2)/1000);
                    tmp2=cat(2,tmp2,tmpHead(1)+tmpHead(2)/1000);
                end
                cellFrame{tmpHead(1)}{tmpHead(2)}.child=tmp;
                cellFrame{tmpHead(1)}{tmpHead(2)}.outflow = cellFrame{tmpHead(1)}{tmpHead(2)}.outflow + xval(i);
                cellFrame{tmpHead(1)}{tmpHead(2)}.cumFlow = cat(2,...
                    cellFrame{tmpHead(1)}{tmpHead(2)}.cumFlow, xval(i));
                
                cellFrame{tmpTail(1)}{tmpTail(2)}.parent = tmp2;
                cellFrame{tmpTail(1)}{tmpTail(2)}.inflow = cellFrame{tmpTail(1)}{tmpTail(2)}.inflow + xval(i);
            end
        end
    end
end