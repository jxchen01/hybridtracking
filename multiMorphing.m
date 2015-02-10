function P = multiMorphing(Ps,Pt,sid)

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

src_idx = find(any(flowMat>1e-5,2));

numMatch = targetLength;
tar_idx = zeros(1,numMatch);

for i=head_idx:1:tail_idx
    tar_idx(i-head_idx+1) = find(flowMat(src_idx(i),:)>1e-5);
end

if(tar_idx(1)>tar_idx(end))
    flag=-1;
elseif(tar_idx(1)<tar_idx(end))
    flag=1;
else
    error('bad optimization');
end

% remove invalid match (not in order)
pts=P2(tar_idx(1),:);
for i=2:1:numMatch-1
    if(sign(tar_idx(i)-tar_idx(i-1))==flag && sign(tar_idx(end)-tar_idx(i))==flag)
        pts=cat(1,pts,P2(tar_idx(i),:));
    end
end
pts=cat(1,pts,P2(tar_idx(end),:));
pts=round(pts);

sz=max([P1(:);P2(:)])+1;
ctl=zeros(sz,sz);
for pid=1:1:size(pts,1)-1
    [px,py]=bresenham(pts(pid,1),pts(pid,2),pts(pid+1,1),pts(pid+1,2));
    pidx = sub2ind([sz,sz],px,py);
    ctl(pidx)=1;
end
ctl = bwmorph(ctl,'thin',Inf);
P=sortOneCellPixel(ctl);