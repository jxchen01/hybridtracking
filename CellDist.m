function contourDiff = CellDist(PsA)

%%%%%%%%%%%
% PsA is the data structure used in contour evolution
%%%%%%%%%%%
nPoints=size(PsA.pts,1);
%%%% fetch PsB %%%%
P=PsA.LastFramePts;
dis=[0;cumsum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)))];
K=zeros(nPoints,2);
K(:,1) = interp1(dis,P(:,1),linspace(0,dis(end),nPoints));
K(:,2) = interp1(dis,P(:,2),linspace(0,dis(end),nPoints));

dv=PsA.pts-K;

na=PsA.normvec;
nb=GetContourNormals2D(K);

d1=abs(dot(dv,na,2));
d2=abs(dot(dv,nb,2));

contourDiff = (mean(d1) + mean(d2))/2;

%%%% visual check %%%%%
% s1=PsA.region; s2=PsB.region;
% tmp=s1+s2*2+1;
% figure,imshow(tmp,rand(4,3));

