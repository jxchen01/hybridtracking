function dist = morphDist(P1, P2, flowCap)

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

[~,fval,exitflag,~] = intlinprog(f,1:length(f),A,b,ones(1,N),...
    flowCap,zeros(1,N),ones(1,N),options);

if(exitflag~=1)
    disp(output.message);
    keyboard
end

dist = fval/flowCap;