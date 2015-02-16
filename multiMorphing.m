function P = multiMorphing(Ps,Pt,sid,sz)

P1=[]; targetLength=0;
kpoint=[];
for i=1:1:numel(Ps)
    kpoint=cat(1,kpoint,size(P1,1)+1);
    pts=Ps{i}.ctl;
    if(i==sid)
        head_idx=size(P1,1)+1;
        targetLength = size(pts,1);
        tail_idx=head_idx+targetLength-1;
    end
    P1=cat(1,P1,pts);
    kpoint=cat(1,kpoint,size(P1,1));
end

if(~targetLength)
    disp('fatal error');
    keyboard
end

P2=Pt.ctl;
flowCap = Pt.inflow;

options = optimoptions('intlinprog','Display', 'off');

len1=size(P1,1);
len2=size(P2,1);

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

for i=1:1:len1
    idx = sub2ind([len1,len2],repmat(i,1,len2),1:1:len2);
    A(i,idx)=1;
end

for i=1:1:len2
    idx = sub2ind([len1,len2],1:1:len1,repmat(i,1,len1));
    A(i+len1, idx)=1;
end

beq=[flowCap;2.0];
Aeq=cat(1,ones(1,N),zeros(1,N)); 
idx=sub2ind([len1,len2],kpoint,ones(numel(kpoint),1));
Aeq(2,idx)=1;
idx=sub2ind([len1,len2],kpoint,repmat(len2,numel(kpoint),1));
Aeq(2,idx)=1;


[xval,~,exitflag,~] = intlinprog(f,1:N,A,ones(M,1),Aeq,beq,zeros(1,N),ones(1,N),options);

if(exitflag~=1)
    disp(output.message);
    keyboard
end

flowMat=zeros(len1,len2);
flowMat(xval>0.5)=1;

tar_idx = [];

for i=head_idx:1:tail_idx
    tt=find(flowMat(i,:)>1e-5);
    if(~isempty(tt))
        tar_idx=cat(1,tar_idx,tt);
    end
end


[m1,midx1] = LIS(tar_idx(1:1:end));
[m2,midx2] = LIS(tar_idx(end:-1:1));

if(m1<m2)
    pts=P2(midx2,:);
else
    pts=P2(midx1,:);
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