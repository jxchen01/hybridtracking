function [m,mlist]=LIS(L)

len = length(L);

lis=zeros(1,len);
lastnode=zeros(1,len);

lis(1)=1;
lastnode(1)=1;
for i=2:1:len
    [t1,lastnode(i)]= maxL(i,L,lis);
    lis(i)=t1+1;
end

[m,idx]=max(lis);
mlist = zeros(1,m);
sig=0;
while(sig<m)
    mlist(end-sig)=L(idx);
    idx=lastnode(idx);
    sig=sig+1;
end

end

function [k, ln] = maxL(idx, a, lis)
    ln=idx;
    for i=1:1:idx-1
        if(a(i)<a(idx) && lis(i)>lis(idx))
            lis(idx) = lis(i);
            ln=i;
        end
    end
    k = lis(idx);
end

