%%%%%%%%%%%% main entry for hybrid tracking %%%%%%%%%%%%%
%%% Created by Jianxu Chen (University of Notre Dame) %%%
%%%%%%%%%%%%%%%%%%% Date: Jan. 2015 %%%%%%%%%%%%%%%%%%%%%
clc
disp('Program Starts...');

sq=5;
RawType='.png';

%%%%% load segmentation results %%%%%%
S=load(['/Users/JianxuChen/Desktop/Research/Myxo_Bacteria/MICCAI2015/data/sq',num2str(sq),'/seg.mat']);
cellEachFrame = S.cellEachFrame;
matEachFrame = S.matEachFrame;

BW = im2bw(imread(['/Users/JianxuChen/Desktop/Research/Myxo_Bacteria/MICCAI2015/data/sq'...
    ,num2str(sq),'/manual.png']));

%%%%% parameters %%%%%
[xdim,ydim]=size(BW);
Options=setParameters(xdim,ydim);
numFrameAhead = Options.numFrameAhead;
numFrame = length(cellEachFrame);
%numFrame = 6;

% load manual segmentation of first frame
ctlImg=bwmorph(BW,'thin',Inf);
cc=bwconncomp(BW);
bwLabel=labelmatrix(cc);
P=ExtractCells(bwLabel,ctlImg,Options);
cellEachFrame{1}=P;
matEachFrame{1}.Mat = double(labelmatrix(bwconncomp(ctlImg)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% low-level association (frame-by-frame matching) %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('low-level starts ...');
srcCellList=cellEachFrame{1};
srcMat=matEachFrame{1}.Mat;
for i=2:1:numFrame
    disp(['Frame: ',num2str(i)])
    tarCellList=cellEachFrame{i};
    tarMat=matEachFrame{i}.Mat;
    [srcCellList,tarCellList]=local_EMD(srcCellList,tarCellList, srcMat, tarMat,Options);
    
    cellEachFrame{i-1}=srcCellList;
    if(i==numFrame)
        cellEachFrame{i}=tarCellList;
    end
    srcCellList=tarCellList;
    srcMat=tarMat; 
end

%%%% main loop %%%
I1=mat2gray(imread(['/Users/JianxuChen/Desktop/Research/Myxo_Bacteria/MICCAI2015/data/sq',...
    num2str(sq),'/raw/img0101',RawType]));
I2=mat2gray(imread(['/Users/JianxuChen/Desktop/Research/Myxo_Bacteria/MICCAI2015/data/sq',...
    num2str(sq),'/raw/img0102',RawType]));

for frameIdx = 2:1:numFrame-numFrameAhead
    disp(['processing frame: ',num2str(frameIdx)]);
    % build correspondence within a period of time
    cellSemiGlobal = Global_EMD(cellEachFrame(1,frameIdx-1:1:frameIdx+numFrameAhead),...
        matEachFrame(1,frameIdx-1:1:frameIdx+numFrameAhead), Options);
    
    % extract (1) confirmed segmentation; (2) confirmed entering cell; 
    % (3) contours needs to evolve
    [Ps,newCellFrame,BMap,propagateIdx] = ConvertCellToContour(cellSemiGlobal,[xdim,ydim],Options);
    
    I3 = mat2gray(imread(['/Users/JianxuChen/Desktop/Research/Myxo_Bacteria/MICCAI2015/data/sq',...
            num2str(sq),'/raw/img0',num2str(100+frameIdx+1),RawType]));
        
    I = mat2gray((I1+I2+I3)./3);

    if(numel(Ps)>0)
        % contour evolution
        [newPs, divisionIDX]=OpenActiveContour(I,Ps,BMap,Options);
    else
        newPs=[]; divisionIDX=[];
    end
    
    % update
    [cellFrame, cMat]=updateCellEachFrame(cellEachFrame(1,frameIdx-1:1:frameIdx+numFrameAhead)...
    ,newCellFrame, newPs, propagateIdx, matEachFrame{frameIdx+1}.Mat ,[xdim,ydim], divisionIDX, Options);

    DrawSegmentedArea2D(cellFrame{2},mat2gray(I),2);
    
    saveas(gcf,['./track/img0',num2str(frameIdx+100),'.png'],'png');
    
    matEachFrame{frameIdx}.Mat = cMat;
    cellEachFrame(1,frameIdx-1:1:frameIdx+1)=cellFrame(1,1:3);
    
    I1=I2;
    I2=I3;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% high-level association (globally min-cost flow) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('high-level starts');
%cellEachFrameGlobal = Global_EMD(cellEachFrame, matEachFrame);
%cellEachFrameLocal = cellEachFrame;
