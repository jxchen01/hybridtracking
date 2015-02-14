function [srcCellList,tarCellList]=local_EMD(srcCellList, tarCellList, srcMat, tarMat, algOptions)

%%%%%%% parameters %%%%%%%
%halfROIs = algOptions.halfROIs;
BoundThreshold = algOptions.BoundThresh;
candiRadius=algOptions.candiRadius;
bodyRatio=algOptions.bodyRatio;
%options = optimset('Algorithm','simplex','Display', 'off', 'Diagnostics',...
%    'off','LargeScale', 'off', 'Simplex', 'on');
options = optimset('Display', 'off', 'Diagnostics','off');

% initialization
srcNum = length(srcCellList);
tarNum = length(tarCellList);

imgSize = size(srcMat);
dimx=imgSize(1);dimy=imgSize(2);

if(size(srcMat)~=size(tarMat))
    error('error in local EMD');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% compute the distance between signitures %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distMat=zeros(srcNum+2,tarNum+2);
validMat=zeros(srcNum,tarNum);

% compute the feasible pairs
tmpSingle=zeros(dimx,dimy);
se = strel('disk', candiRadius);
for i=1:1:srcNum
    cellA=srcCellList{i}.ctl;
    tmpIdx = sub2ind(imgSize,cellA(:,1),cellA(:,2));
    tmpSingle(tmpIdx)=1;

    % find overlaping regions
    tmpSingleRegion=imdilate(tmpSingle,se);
    tmpCandi = unique(nonzeros(tmpSingleRegion.*tarMat));
    validMat(i,tmpCandi)=1;
    
    % reset tmpSingle
    tmpSingle(tmpIdx)=0;
end
clear tmpSingle se tmpCandi tmpIdx cellA

% compute the distance between each feasible pair
for i=1:1:srcNum
    sid = i;
    tmpValid = find(validMat(i,:)>0);
    for j=1:1:length(tmpValid)
        tid = tmpValid(j);
        tmpD = ComputeDist(srcCellList{sid}.ctl,srcCellList{sid}.length,...
                           tarCellList{tid}.ctl,tarCellList{tid}.length,0);       
        if(tmpD<500)
            distMat(sid,tid)=tmpD;
            srcCellList{sid}.candi=cat(1,srcCellList{sid}.candi,tid);
        else
            validMat(sid,tid)=0;
        end
    end
end
clear tid sid tmpD tmpValid 

%%%%%%%%%%%% compute the possibility of leaving %%%%%%%%%%%%%%%
% stort in distMat(:, tarNum+1)
for i=1:1:srcNum
    lenA = srcCellList{i}.length;
    cellA = srcCellList{i}.ctl;
    
    xmin=min([cellA(1,1), dimx-cellA(1,1), cellA(end,1), dimx-cellA(end,1)]);
    ymin=min([cellA(1,2), dimy-cellA(1,2), cellA(end,2), dimy-cellA(end,2)]);
    dmin=min([xmin,ymin]);
    
    if(dmin<BoundThreshold)
        distMat(i,tarNum+1)=0.43*(dmin+lenA/4); % 30-degree
        srcCellList{i}.relaxOutCost=distMat(i,tarNum+1);
    else
        %distMat(i,tarNum+2)=0.75*0.5*srcCellList{i}.length;
        distMat(i,tarNum+2)=0.8*40; % assume length=80
        srcCellList{i}.relaxOutCost=-1*distMat(i,tarNum+2);
    end
end

%%%%%%%%%%%%%%% compute the possibility of entering %%%%%%%%%%%%%%%
% store in distMat(srcNum+1, :)
for i=1:1:tarNum
    lenB=tarCellList{i}.length;
    cellB=tarCellList{i}.ctl;

    xmin=min([cellB(1,1), dimx-cellB(1,1), cellB(end,1), dimx-cellB(end,1)]);
    ymin=min([cellB(1,2), dimy-cellB(1,2), cellB(end,2), dimy-cellB(end,2)]);
    dmin=min([xmin,ymin]);
    
    if(dmin<BoundThreshold)
        distMat(srcNum+1,i)=0.43*(dmin+lenB/4); % 30-degree
        tarCellList{i}.relaxInCost=distMat(srcNum+1,i);
    else
        %distMat(srcNum+2,i)=0.75*0.5*tarCellList{i}.length;
        distMat(srcNum+2,i)=0.8*40; %assume length=80
        tarCellList{i}.relaxInCost=-1*distMat(srcNum+2,i);
    end
end

clear lenA lenB dmin hmin tmin xmin ymin cellA cellB lenA lenB i j im1 im2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% perform matching on the whole region %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sList=1:1:srcNum;
sNum=srcNum;

% get corresponding tL
tmp=zeros(1,tarNum);
varNum=0;
for i=1:1:sNum
    varNum = varNum + nnz(validMat(sList(i),:));
    tmp=tmp+validMat(sList(i),:);
end
tList=find(tmp>0);
tNum=numel(tList);
clear tmp

% get leaving variables and relaxation variables
varRS=zeros(sNum,1);
for i=1:1:sNum
    if(distMat(sList(i),tarNum+1)>0)
        varRS(i)=distMat(sList(i),tarNum+1);
    else
        varRS(i)=distMat(sList(i),tarNum+2);
    end
end

% get entering variables and relaxation variables
varRT=zeros(tNum,1);
for i=1:1:tNum
    if(distMat(srcNum+1,tList(i))>0)
        varRT(i)=distMat(srcNum+1,tList(i));
    else
        varRT(i)=distMat(srcNum+2,tList(i));
    end
end

% build the linear optimization problem
weight=zeros(1,varNum+sNum+tNum);
A = zeros(sNum+tNum,varNum+sNum+tNum);  % flow constraint at each node
b = zeros(sNum+tNum,1);
lb = zeros(varNum+sNum+tNum,1);
ub = [];
S1=0;
sig=0;
% update b, S1, weight and A w.r.t. srcMat
for i=1:1:sNum
    b(i)=srcCellList{i}.length;
    S1=S1+b(i);
    
    tmpNbr = find(validMat(i,:)>0);
    for j=1:1:length(tmpNbr)
        sig=sig+1;
        weight(sig)=distMat(i,tmpNbr(j));
        A(i,sig)=1;
        tListInd=find(tList==tmpNbr(j));
        A(tListInd+sNum,sig)=1;
    end
    A(i,varNum+i)=1;
end
clear tListInd

weight(1,sig+1:sig+sNum)=varRS(:,1);
sig=sig+sNum;
weight(1,sig+1:sig+tNum)=varRT(:,1);

S2=0;
for i=1:1:tNum
    b(i+sNum)=tarCellList{tList(i)}.length;
    S2=S2+b(i+sNum);
    A(i+sNum,varNum+sNum+i)=1;
end

% build Aeq, beq
Aeq = ones(1,varNum+sNum+tNum);
beq = min([S1,S2]);

[xval,~,exitflag,output] = linprog(weight,A,b,Aeq,beq,lb,ub,[], options);
if(exitflag~=1)
    disp(output.message);
    keyboard
end

%%% feasible matching
matchingMat=zeros(srcNum,tarNum);
sig=0;
for i=1:1:sNum
    tmpNbr = find(validMat(i,:)>0);
    for j=1:1:numel(tmpNbr)
        sig=sig+1;
        tf=xval(sig);
        if(tf>0) % posive flow !
            indInTar=tmpNbr(j);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Do local matching conservatively %%%%%%%
        %%%%% i.e. only keep very good matching %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(distMat(i,indInTar)>3)
            continue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % treat positve flow as a matching, only if
        % the flow amount if greater than the length
        % of either src or tar
        if(tf>bodyRatio*srcCellList{i}.length || tf>bodyRatio*tarCellList{indInTar}.length)
            matchingMat(i,indInTar)=tf;
        end
    end
end

clear indInTar tf tmpNbr sig len

% update parent, child, inflow and outflow
for i=1:1:sNum
    sid=sList(i);
    new_child=find(matchingMat(sid,:)>0);
    if(numel(new_child)==0)
        continue;
    end
    
    srcCellList{sid}.child=new_child;
    for j=1:1:numel(new_child)
        tmpFlow = min([tarCellList{new_child(j)}.length-tarCellList{new_child(j)}.inflow,srcCellList{sid}.length-srcCellList{sid}.outflow]);
        srcCellList{sid}.outflow = srcCellList{sid}.outflow + tmpFlow;
        srcCellList{sid}.cumFlow = cat(2,srcCellList{sid}.cumFlow, tmpFlow);
        
%         if(isempty(tarCellList{new_child(j)}.parent))
%             tarCellList{new_child(j)}.parent = sid;
%             tarCellList{new_child(j)}.inflow = tmpFlow;
%         else
            tarCellList{new_child(j)}.parent=...
                cat(2,tarCellList{new_child(j)}.parent,sid);
            tarCellList{new_child(j)}.inflow = tarCellList{new_child(j)}.inflow + tmpFlow;
%         end      
    end
        
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% perform matching on each ROI %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROImask=zeros(dimx,dimy);
% countRegion=0;
% for rx1=1:halfROIs:dimx
%     [rx2,flagx]=min([2*halfROIs+rx1,dimx]);
%     for ry1=1:halfROIs:dimy    
%         countRegion=countRegion+1;
%         % disp(['Region: ',num2str(countRegion)]);       
%         [ry2,flagy]=min([2*halfROIs+ry1,dimy]);
%         ROImask(rx1:rx2,ry1:ry2)=1;
%         tmp=srcMat.*ROImask;
%         sList=unique(nonzeros(tmp));
%         sNum=numel(sList);
%         clear tmp
%         
%         % get corresponding tL 
%         tmp=zeros(1,tarNum);
%         varNum=0;
%         for i=1:1:sNum
%             varNum = varNum + nnz(validMat(sList(i),:));
%             tmp=tmp+validMat(sList(i),:);
%         end
%         tList=find(tmp>0);
%         tNum=numel(tList);
%         clear tmp
%         
%         % get leaving variables and relaxation variables
%         varRS=zeros(sNum,1);
%         for i=1:1:sNum
%             if(distMat(sList(i),tarNum+1)>0)
%                 varRS(i)=distMat(sList(i),tarNum+1);
%             else
%                 varRS(i)=distMat(sList(i),tarNum+2);
%             end
%         end   
%         
%         % get entering variables and relaxation variables
%         varRT=zeros(tNum,1);
%         for i=1:1:tNum
%             if(distMat(srcNum+1,tList(i))>0)
%                 varRT(i)=distMat(srcNum+1,tList(i));
%             else
%                 varRT(i)=distMat(srcNum+2,tList(i));
%             end
%         end
% 
%         % build the linear optimization problem
%         weight=zeros(1,varNum+sNum+tNum);
%         A = zeros(sNum+tNum,varNum+sNum+tNum);  % flow constraint at each node
%         b = zeros(sNum+tNum,1);
%         lb = zeros(varNum+sNum+tNum,1);
%         ub = [];
%         S1=0;
%         sig=0;
%         for i=1:1:sNum
%             indInSrc=sList(i);
%             tmpNbr = find(validMat(indInSrc,:)>0);
%             len = length(tmpNbr);
%             
%             b(i)=srcCellList{indInSrc}.length;
%             S1=S1+b(i);
%             
%             for j=1:1:len
%                 sig=sig+1;
%                 indInTar=tmpNbr(j);
%                 weight(sig)=distMat(indInSrc,indInTar);
%                 A(i,sig)=1;
%                 tListInd=find(tList==indInTar);
%                 A(tListInd+sNum,sig)=1;
%             end
%             A(i,varNum+i)=1;
%         end
%         clear tListInd
% 
%         weight(1,sig+1:sig+sNum)=varRS(:,1);
%         sig=sig+sNum;
%         weight(1,sig+1:sig+tNum)=varRT(:,1);
%         
%         S2=0;
%         for i=1:1:tNum
%             indInTar=tList(i);
%             b(i+sNum)=tarCellList{indInTar}.length;
%             S2=S2+b(i+sNum);
%             A(i+sNum,varNum+sNum+i)=1;
%         end
%         Aeq = ones(1,varNum+sNum+tNum);
%         beq = min([S1,S2]);
%         
%         [xval,~,exitflag,output] = linprog(weight,A,b,Aeq,beq,lb,ub,[], options);
%         if(exitflag~=1)
%             disp(output.message);
%             keyboard
%         end
%         
%         %%% feasible matching
%         matchingMat=zeros(srcNum,tarNum);
%         sig=0;
%         for i=1:1:sNum
%             indInSrc=sList(i);
%             tmpNbr = find(validMat(indInSrc,:)>0);
%             len = length(tmpNbr);
%             for j=1:1:len
%                 sig=sig+1;
%                 tf=xval(sig);
%                 if(tf>0)
%                     indInTar=tmpNbr(j);
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%% What ??? %%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 if(distMat(indInSrc,indInTar)>1)
% %                     continue;
% %                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 if(tf>bodyRatio*srcCellList{indInSrc}.length || tf>bodyRatio*tarCellList{indInTar}.length)
%                     matchingMat(indInSrc,indInTar)=tf;
%                 end
%             end
%         end    
% %             % check sList
% %             if(numel(find(matchingMat(indInSrc,:)>0))>1)
% %                 tt=find(matchingMat(indInSrc,:)>srcCellList{indInSrc}.length);
% %                 if(numel(tt)>0)
% %                     keyboard % why? for debug and check
% %                     tp=matchingMat(indInSrc,tt);
% %                     matchingMat(indInSrc,:)=0;
% %                     matchingMat(indInSrc,tt)=tp;
% %                 end
% %             end
% %         % check tList
% %         for i=1:1:tNum
% %             indInTar=tList(i);
% %             if(numel(find(matchingMat(:,indInTar)>0))>1)
% %                 tt=find(matchingMat(:,indInTar)>tarCellList{indInTar}.length);
% %                 if(numel(tt)>0)
% %                     keyboard % why? for debug and check
% %                     tp=matchingMat(tt,indInTar);
% %                     matchingMat(:,indInSrc)=0;
% %                     matchingMat(tt,indInTar)=tp;
% %                 end
% %             end
% %         end
%         clear indInSrc indInTar tf tmpNbr sig len
%         
%         %%% compare with existing matching results (from other ROI)
%         for i=1:1:sNum
%             sid=sList(i);
%             new_child=find(matchingMat(sid,:)>0);
%             if(numel(new_child)==0)
%                 continue;
%             end
%             len=numel(srcCellList{sid}.child);           
%             if(len>0)
%                 flag=0;
%                 if(len==length(new_child))
%                     diffNum=numel(setxor(srcCellList{sid}.child,new_child));
%                     if(diffNum>0)
%                         flag=1;
%                     end
%                 else
%                     flag=1;
%                 end
%                 if(flag>0)
%                     %%%% for debug
%                     % keyboard
%                     %%%% 
%                     trm=srcCellList{sid}.child;
%                     for j=1:1:len  % len=numel(srcCellList{sid}.child);
%                        
%                         tid=trm(j);
%                         if(tid==0)
%                             continue;
%                         end
%                         %%% for debug
%                         %if(numel(tarCellList{tid}.parent)>1)
%                         %    keyboard;
%                         %end
%                         %%% 
%                         tmpParent=tarCellList{tid}.parent;
%                         tarCellList{tid}.parent=tmpParent(tmpParent~=sid);
%                         tarCellList{tid}.inflow = tarCellList{tid}.inflow-flowMat(sid,tid); 
%                         flowMat(sid,tid)=0;
%                     end
%                     srcCellList{sid}.child=0;
%                     srcCellList{sid}.outflow=0;
%                 end
%             else
%                 srcCellList{sid}.child=new_child;
%                 for j=1:1:numel(new_child)
%                     tmpFlow = min([tarCellList{new_child(j)}.length-tarCellList{new_child(j)}.inflow,srcCellList{sid}.length-srcCellList{sid}.outflow]);
%                     srcCellList{sid}.outflow = srcCellList{sid}.outflow + tmpFlow;
%                         
%                     if(isempty(tarCellList{new_child(j)}.parent))
%                         tarCellList{new_child(j)}.parent = sid;
%                         tarCellList{new_child(j)}.inflow = tmpFlow;
%                     else
%                         tarCellList{new_child(j)}.parent=...
%                             cat(2,tarCellList{new_child(j)}.parent,sid);
%                         tarCellList{new_child(j)}.inflow = tarCellList{new_child(j)}.inflow + tmpFlow;
%                     end
%                     
%                     try
%                         flowMat(sid,new_child(j))=tmpFlow;
%                     catch
%                         keyboard
%                     end
%                 end
%             end
%                 
%         end
%         
%         if(flagy==2)
%             break;  % stop moving ROI in y direction
%         end
%     end 
%     if(flagx==2)
%         break;      % stop moving ROI in x direction
%     end
% end
% 
% clear fval xval i j lb ub Aeq beq weight A b options new_child 
% clear sid tid tt S1 S2 diffNum exitflag flag flagx flagy tmpFlow