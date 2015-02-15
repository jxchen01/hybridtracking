function P=cellMorphing(P1,P2,flowCap,ftime,sz)

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

if(ftime==1)
    pts=P2(midxb,:);
else
    pts=zeros(mb,2);
    for i=1:1:mb
        sid = find(tar_idx==midxb(i));
        if(isempty(sid))
            disp('error in EMD morphing');
            keyboard;
        end
        pts(i,:) = P1(src_idx(sid),:) + (P2(midxb(i),:) - P1(src_idx(sid),:))./double(ftime);
    end
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
P=sortOneCellPixel(ctl);


