function [cellList, cMat] = ExtractCells(bwLabel,im,Options)

im(im>0)=1;
[xdim,ydim]=size(im);
nbr=[1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];

epImg = bwmorph(im, 'endpoint');
labelImg=epImg+im;
cellNum = nnz(epImg)/2;
if(mod(nnz(epImg),2)==1)
    disp('error in number of endpoints');
    error('error in extracting manual seg');
end
clear epImg

ep=find(labelImg==2);
currentCellNum=0;
cellList=cell(1,cellNum);
cMat=zeros(xdim,ydim);

while(~isempty(ep))
    % extract a new end point
    [xt,yt]=ind2sub([xdim ydim],ep(1));
    labelImg(xt,yt)=0;
    
    oep=0; % flag indicating if other endpoint is removed in iteration
    
    % extract the cell body
    tmpPixInd=zeros(50,2);
    pixNum=1;
    tmpPixInd(1,1)=xt; tmpPixInd(1,2)=yt;
    flag=1;
    while(flag)
        flag=0;
        for i=1:8
            x0=xt+nbr(i,1);
            y0=yt+nbr(i,2);
            if(x0<1 || x0>xdim || y0<1 || y0>ydim)
                continue;
            end
            if(labelImg(x0,y0)>0)
                flag=1;
                if(labelImg(x0,y0)==2)
                    if(oep==0)
                        oep=1;
                    else
                        disp('more then two endpoints are removed in one iteration');
                        error('error in extracting manual seg');
                    end
                end
                
                pixNum=pixNum+1;
                tmpPixInd(pixNum,1)=x0;tmpPixInd(pixNum,2)=y0;
                labelImg(x0,y0)=0;
                xt=x0; yt=y0;
                
                break;
            end
        end
    end
    
    if(oep==0)
        disp('only one endpoint is removed in this iteration');
        error('error in extracting manual seg');
    end
    
    %%%%% remove short cells %%%%
    if(pixNum<Options.pixelCanSkip)
        ep=find(labelImg==2);
        continue;
    end
    
    currentCellNum=currentCellNum+1;   
    pts=zeros(pixNum,2);
    pts(:,:)=tmpPixInd(1:pixNum,:);    
    seg_region = (bwLabel==bwLabel(pts(1,1),pts(1,2)));
    
    cellList{currentCellNum}=struct('length',pixNum,'ctl',pts,'child',[],...
        'parent',[],'candi',[],'inflow',0,'outflow',0,'relaxInCost',0,...
        'relaxOutCost',0,'seg',seg_region,'id',currentCellNum,'cumFlow',[]);
    
    if(~isCloseToBoundary(pts,xdim,ydim, Options.BoundThresh))
        cellList{currentCellNum}.copyLength = pixNum;
    end

    pts_idx = sub2ind([xdim,ydim],pts(:,1),pts(:,2));
    cMat(pts_idx)=currentCellNum;
    
    ep=find(labelImg==2); 
end