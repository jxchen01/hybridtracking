function im1=pruneLine(im1)

%%%%% check if pruning is necessary
bp=bwmorph(im1,'branchpoints');
if(nnz(bp)==0)
    return;
end

%%%%% start pruning
MinBranch=5;

CopyImg=im1;  % make a copy for debug
im1(im1>0)=1;
[xdim,ydim]=size(im1);

mask=bwmorph(bp,'dilate');
mask(bp>0)=0;
mask(im1==0)=0; % contains the points connecting to branchpoints

% mark branch point as 20 and points next to bp (on branches) as 10
labelImg=im1+9*mask;
labelImg=labelImg+19*bp;

nbr=[1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];

cp=find(labelImg==10); % those points next to bp

while(~isempty(cp))
    [xt,yt]=ind2sub([xdim ydim],cp(1));
    labelImg(xt,yt)=0;
    
    tmpPixInd=zeros(MinBranch,2);
    pixNum=1;
    tmpPixInd(1,1)=xt; tmpPixInd(1,2)=yt;
    xp=xt; yp=yt; % previous postion
    
    flag=1;
    while(flag>0)
        flag=0;
        for i=1:8
            x0=xt+nbr(i,1);
            y0=yt+nbr(i,2);
            if(x0<1 || x0>xdim || y0<1 || y0>ydim)
                continue;
            end
            if(labelImg(x0,y0)==1 && (x0~=xp || y0~=yp) )
                pixNum=pixNum+1;
                tmpPixInd(pixNum,1)=x0;tmpPixInd(pixNum,2)=y0;
                xp=xt; yp=yt;
                xt=x0; yt=y0;
                if(pixNum==MinBranch)
                    flag=-1;
                else
                    flag=1;
                end
                break;
            end
        end
    end
    
    if(flag==0) % this is a branch shorter than MinBranch
        for i=1:1:pixNum
            im1(tmpPixInd(i,1),tmpPixInd(i,2))=0;
        end    
        % remove the branchpoint when the remaining branches keep connected
        [xt,yt]=ind2sub([xdim ydim],cp(1));
        flag=0;
        for i=1:1:8
            x0=xt+nbr(i,1);
            y0=yt+nbr(i,2);
            if(x0<1 || x0>xdim || y0<1 || y0>ydim)
                continue;
            end
            if(bp(x0,y0)>0)
                flag=1;
                break;
            end
        end
        if(flag>0)
            ocp=zeros(3,2);
            ocpNum=0;
            for i=1:1:8
                x1=x0+nbr(i,1);
                y1=y0+nbr(i,2);
                if(x1<1 || x1>xdim || y1<1 || y1>ydim)
                    continue;
                end
                if(labelImg(x1,y1)==10)
                    ocpNum=ocpNum+1;
                    ocp(ocpNum,1)=x1;
                    ocp(ocpNum,2)=y1;
                end
            end
            if(ocpNum>2)
                disp('degree > 3. Recheck the code');
                keyboard
            elseif(ocpNum==2)
                if(abs(ocp(1,1)-ocp(2,1))<2 && abs(ocp(1,2)-ocp(2,2))<2)
                    im1(x0,y0)=0;  %%%% remove the branch point
                end
            end
        else
            disp('error in removing bp');
        end
    end
    
    cp=find(labelImg==10);
end

cc=bwconncomp(im1);
if(cc.NumObjects~=1)
    disp('cell is broken');
    keyboard
end

end
