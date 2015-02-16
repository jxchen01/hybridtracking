function P=cellMorphing(P1,P2,flowCap,ftime,sz)

if(size(P1,2)~=2 || size(P2,2)~=2)
    disp('error in dimenstion');
    keyboard;
end

len1=size(P1,1);
len2=size(P2,1);

if(ftime==1)
    options = optimoptions('intlinprog','Display', 'off');
    distMat = zeros(len1,len2);
    for i=1:1:len1
        for j=1:1:len2
            distMat(i,j)=norm(P1(i,:)-P2(j,:));
        end
    end

    f=distMat(:);
    N=numel(f);
    M=len1+len2;
    A=zeros(M,N);
    b=ones(M,1);

    for i=1:1:len1
        idx = sub2ind([len1,len2],repmat(i,1,len2),1:1:len2);
        A(i,idx)=1;
    end
    
    for i=1:1:len2
        idx = sub2ind([len1,len2],1:1:len1,repmat(i,1,len1));
        A(i+len1, idx)=1;
    end
    
    [xval,fval,exitflag,~] = intlinprog(f,1:length(f),A,b,ones(1,N),...
        double(flowCap),zeros(1,N),ones(1,N),options);

    if(exitflag~=1)
        disp(output.message);
        error('error in cellMorhing optimization');
    end

    flowMat=zeros(len1,len2);
    flowMat(xval>0.5)=1;
    
    src_idx = [];
    tar_idx = [];
    
    for i=1:1:len1
        tt=find(flowMat(i,:)>1e-5);
        if(~isempty(tt))
            tar_idx=cat(1,tar_idx,tt);
            src_idx=cat(1,src_idx,i);
        end
    end
    
    [m1,midx1] = LIS(tar_idx(1:1:end));
    [m2,midx2] = LIS(tar_idx(end:-1:1));
    
    if(m1<m2)
        mb=m2; midxb=midx2;
    else
        mb=m1; midxb=midx1;
    end
    
    pts=P2(midxb,:);

    %     pts=zeros(mb,2);
    %     for i=1:1:mb
    %         sid = find(tar_idx==midxb(i));
    %         if(isempty(sid))
    %             disp('error in EMD morphing');
    %             keyboard;
    %         end
    %         pts(i,:) = P1(src_idx(sid),:) + (P2(midxb(i),:) - P1(src_idx(sid),:))./double(ftime);
    %     end
    
    

elseif(ftime>1.5)
    if(len1<=len2)
        LA=len1; LB=len2;
        A=P1; B=P2;
    else
        LA=len2; LB=len1;
        A=P2; B=P1;
    end
    
    [MaxD1, best_p1]=cellAlign(A, LA, B, LB);
    [MaxD2, best_p2]=cellAlign(A, LA, fliplr(B), LB);
    
    if(MaxD1>MaxD2)
        B=fliplr(B);
        k=best_p2;
    else
        k=best_p1;
    end
    
    pts=zeros(LA,2);
    pts(:,:) = A(1:LA,:) + (B((1+k):(LA+k),:)-A(1:LA,:))./double(ftime);

    
else
    error('error in ftime');
end

pts=round(pts);
pts(pts(:,1)<1,:)=1; pts(pts(:,1)>sz(1),:)=sz(1);
pts(pts(:,2)<1,:)=1; pts(pts(:,2)>sz(2),:)=sz(2);

ctl=zeros(sz);
for pid=1:1:size(pts,1)-1
    [px,py]=bresenham(pts(pid,1),pts(pid,2),pts(pid+1,1),pts(pid+1,2));
    pidx = sub2ind(sz,px,py);
    ctl(pidx)=1;
end
ctl = bwmorph(ctl,'thin',Inf);
    
bp_test=bwmorph(ctl,'branchpoint');
if(any(bp_test(:)))
    se0=strel('disk',3,0);
    cr= imdilate(ctl,se0);
    ctl=bwmorph(cr,'thin',Inf);
    clear cr se0
end
clear bp_test

P=sortOneCellPixel(ctl);

end



function [MaxD, flag]=cellAlign(A, LA, B, LB)

    MaxD=Inf;
    flag=0;
    tM=zeros(LA,2);
    td=zeros(LA,1);
    for k=0:1:LB-LA
        tM(1:LA,:)=A(1:LA,:)-B((1+k):(LA+k),:);
        td(1:LA,1)=tM(1:LA,1).*tM(1:LA,1)+tM(1:LA,2).*tM(1:LA,2);
        td=sqrt(td);
        tmpD=sum(td(:));
        tmpD=tmpD/LA;
        if(tmpD<MaxD)
            MaxD=tmpD;
            flag=k;
        end
    end
end
