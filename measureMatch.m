function [MaxD,div]=measureMatch(A,LA, B,LB, debug)
% find the position of best matching
% Assumption: lenA <= lenB
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

% best position is k
tM(1:LA,:)=B((1+flag):(LA+flag),:)-A(1:LA,:); %motion direction
tV=zeros(LA,2); %object1 orientation
tU=zeros(LA,2); %object2 orientation
%tV(1,:)=A(1,:)-A(2,:);
%tV(LA,:)=A(LA-1,:)-A(LA,:);
tV(2:(LA-1),:)=A(1:(LA-2),:)-A(3:LA,:);
tU(2:(LA-1),:)=B((1+flag):(LA+flag-2),:)-B((3+flag):(LA+flag),:);

div=0;
for k=2:1:(LA-1)
    v0=[tM(k,1),tM(k,2)];
    v1=[tV(k,1),tV(k,2)];
    v2=[tU(k,1),tU(k,2)];
    if(norm(v0)==0)
        continue;
    else
        div=div+acos(abs(dot(v1,v0)/(norm(v1)*norm(v0))));
        div=div+acos(abs(dot(v2,v0)/(norm(v2)*norm(v0))));
    end
end
div=cos(0.5*div/(LA-2));
% div=0;
% for k=2:1:(LA-1)
%     v1=[tM(k,1),tM(k,2)];
%     v2=[tV(k,1),tV(k,2)];
%     if(norm(v1)==0)
%         div=div+1;
%     else
%         div=div+abs(dot(v1,v2)/(norm(v1)*norm(v2)));
%     end
% end
% div=div/(LA-2);

if(debug)
    keyboard;
end

