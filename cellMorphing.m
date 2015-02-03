function P=cellMorphing(P1,P2,flowCap,ftime)

if(size(P1,2)~=2 || size(P2,2)~=2)
    disp('error in dimenstion');
    keyboard;
end


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
    keyboard
end

flowMat=zeros(len1,len2);
flowMat(xval>0)=1;

src_idx = find(any(flowMat>1e-5,2));

% remove invalid match (not in order)
numMatch = numel(src_idx);
tar_idx = zeros(1,numMatch);
for i=1:1:numMatch
    tar_idx(i) = find(flowMat(src_idx(i),:)>1e-5);
end

if(tar_idx(1)>tar_idx(end))
    flag=-1;
elseif(tar_idx(1)<tar_idx(end))
    flag=1;
else
    error('bad optimization');
end

pts=P1(src_idx(1),:)+(P2(tar_idx(1),:)-P1(src_idx(1),:))./double(ftime);
for i=2:1:numMatch-1
    if(sign(tar_idx(i)-tar_idx(i-1))==flag && sign(tar_idx(end)-tar_idx(i))==flag)
        pts=cat(1,pts,P1(src_idx(i),:)+(P2(tar_idx(i),:)-P1(src_idx(i),:))./double(ftime));
    end
end
pts=cat(1,pts,P1(src_idx(end),:)+(P2(tar_idx(end),:)-P1(src_idx(end),:))./double(ftime));
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


