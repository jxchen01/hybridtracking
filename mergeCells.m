function P = mergeCells(Ps)

numCell = numel(Ps);
if(numCell==2)
    c1=Ps{1}.ctl;
    c2=Ps{2}.ctl;
    P=combineTwoCells(c1,c2);
else
    % select the closest two cells
    minDist=10000;
    for i=1:1:numCell
        c1=Ps{i}.ctl;
        for j=i+1:1:numCell
            c2=Ps{j}.ctl;
            tD = min([norm(c1(1,:)-c2(1,:))+norm(c1(1,:)-c2(end,:)),...
                norm(c1(end,:)-c2(1,:)), norm(c1(end,:)-c2(end,:))]);
            if(tD<minDist)
                minDist=tD;
                mi=i;mj=j;
            end
        end
    end
    c1=Ps{mi}.ctl;
    c2=Ps{mj}.ctl;   
    P=combineTwoCells(c1,c2);
end

end

function P=combineTwoCells(c1,c2)
    h=zeros(2,2);
    h(1,1)=norm(c1(1,:)-c2(1,:));
    h(1,2)=norm(c1(1,:)-c2(end,:));
    h(2,1)=norm(c1(end,:)-c2(1,:));
    h(2,2)=norm(c1(end,:)-c2(end,:));
    
    [m1, t1]=min(h);
    [~, t2]=min(m1);
    t1=t1(t2);
    
    if(t1==1)
        P=c1(end:-1:1,:);
    else
        P=c1(1:1:end,:);
    end
    
    if(t2==1)
        P=cat(1,P,c2(1:1:end,:));
    else
        P=cat(1,P,c2(end:-1:1,:));
    end
end