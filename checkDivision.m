function divisionIDX = checkDivision(Ps, I)

sz=size(I);
divisionIDX=[];

for i=1:1:numel(Ps)
    P=Ps{i}.pts;
    len=size(P,1);
    x=1:1:len;
    
    t=interp2(I,P(:,2),P(:,1));
    
    y = smooth(x,t,0.2,'loess');

    [vmin,pos]=min(y);
    if(pos>0.4*len && pos<0.6*len && (max(t))/vmin>2.25)
        divisionIDX = cat(2,divisionIDX,i);
    end
            
end
