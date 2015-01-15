function dist=ComputeDist(cellA,lenA, cellB,lenB,debug)
distThresh=30; % caution: not divThresh
%distOffset=0.01;
%coeff=(pi/2-pi/25)/(1-divThresh);

hxA=cellA(1,1); hyA=cellA(1,2);
txA=cellA(lenA,1); tyA=cellA(lenA,2);
hxB=cellB(1,1); hyB=cellB(1,2);
txB=cellB(lenB,1); tyB=cellB(lenB,2);

hh=sqrt((hxA-hxB)*(hxA-hxB)+(hyA-hyB)*(hyA-hyB));
tt=sqrt((txA-txB)*(txA-txB)+(tyA-tyB)*(tyA-tyB));
ht=sqrt((hxA-txB)*(hxA-txB)+(hyA-tyB)*(hyA-tyB));
th=sqrt((txA-hxB)*(txA-hxB)+(tyA-hyB)*(tyA-hyB));

apDist=min([hh,tt,ht,th]);

if(apDist<distThresh)
    tB=zeros(lenB,2);
    for i=1:lenB
        tB(i,:)=cellB(lenB-i+1,:);
    end
    
    if(lenA<lenB)
        if(debug)
            [MaxD1,div1]=measureMatch(cellA,lenA, cellB,lenB,1);
            [MaxD2,div2]=measureMatch(cellA,lenA, tB,lenB,1);
        else
            [MaxD1,div1]=measureMatch(cellA,lenA, cellB,lenB,0);
            [MaxD2,div2]=measureMatch(cellA,lenA, tB,lenB,0);
        end
    else
        if(debug)
            [MaxD1,div1]=measureMatch(cellB,lenB, cellA,lenA,1);
            [MaxD2,div2]=measureMatch(tB,lenB, cellA,lenA,1);
        else
            [MaxD1,div1]=measureMatch(cellB,lenB, cellA,lenA,0);
            [MaxD2,div2]=measureMatch(tB,lenB, cellA,lenA,0);
        end
    end
    
    if(MaxD1<MaxD2) % div=cos(t), so small div means big deviation
%         if(div1>divThresh)
%             dist=MaxD1*cos((div1-divThresh)*coeff);
%         else
%             dist=MaxD1;
%         end
        dist=MaxD1*(1-(div1)^4);
    else
%         if(div2>divThresh)
%             dist=MaxD2*cos((div2-divThresh)*coeff);
%         else
%             dist=MaxD2;
%         end
         dist=MaxD2*(1-(div2)^4);
    end
    if(debug)
        keyboard
    end
else
    dist=10000;
end
