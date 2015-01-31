%%%%%%%%%%%%%%%%%%% main entry for hybrid tracking %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Created by Jianxu Chen (University of Notre Dame) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Date: Jan. 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
disp('Program Starts...');

% data parameters
sq=41;
numFrame=9;
RawType='.png';
fpath = '/Users/JianxuChen/Dropbox/Private/miccai2015/';
% '/Users/JianxuChen/Desktop/Research/Myxo_Bacteria/MICCAI2015/data/'

% load manual segmentation and parameters
BW = im2bw(imread([fpath,'sq',num2str(sq),'/manual.png']));

[xdim,ydim]=size(BW);
Options=setParameters(xdim,ydim);

numFrameAhead = Options.numFrameAhead;
cMap = rand(1000,3).*0.9 + 0.1; % displaying the trackig results
cMap(1,:)=[0,0,0];

BW=regionRefine(BW);
cc=bwconncomp(BW);
bwLabel=labelmatrix(cc);
ctlImg=bwmorph(BW,'thin',Inf);
[P, cMat]=ExtractCells(bwLabel,ctlImg,Options);

maxID = numel(P); % hold the maximum of occupied index value

clear cc bwLabel ctlImg BW 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% initialize the data for main loop %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellBlock=cell(1,numFrameAhead+2);
matBlock =cell(1,numFrameAhead+2);
cellBlock{1}=P;
matBlock{1}=struct('Mat',cMat);
for i=2:1+numFrameAhead
    S=load([fpath,'sq',num2str(sq),'/seg_data/seg0',num2str(100+i),'.mat']);
    matBlock{i}=S.matFrame;
    
    [cellBlock{i-1},cellBlock{i}]=local_EMD(cellBlock{i-1},S.cellFrame,...
        matBlock{i-1}.Mat,matBlock{i}.Mat,Options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%     main loop     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I1=mat2gray(imread([fpath,'sq',num2str(sq),'/raw/img0101',RawType]));
I2=mat2gray(imread([fpath,'sq',num2str(sq),'/raw/img0102',RawType]));

for frameIdx = 2:1:numFrame-numFrameAhead
    disp(['processing frame: ',num2str(frameIdx)]);
    
%     if(frameIdx==14)
%         keyboard;
%     end
    
    % add the new Frame, index = frameIdx + numFrameAhead
    indNew = frameIdx + numFrameAhead;
    S=load([fpath,'sq',num2str(sq),'/seg_data/seg0',num2str(100+indNew),'.mat']);
    matBlock{numFrameAhead+2} = S.matFrame;
    [cellBlock{numFrameAhead+1},cellBlock{numFrameAhead+2}]=...
        local_EMD(cellBlock{numFrameAhead+1},S.cellFrame,...
        matBlock{numFrameAhead+1}.Mat,matBlock{numFrameAhead+2}.Mat,Options);
    
    clear S
    
    % build the image of interest
    I3 = mat2gray(imread([fpath,'sq',num2str(sq),'/raw/img0',num2str(100+frameIdx+1),RawType]));    
    %I = mat2gray((I1+I2+I3)./3);
    I=I2;
    
    % build correspondence within a period of time
    cellSemiGlobal = Global_EMD(cellBlock, matBlock, Options);
    
    % extract:
    %   (1) confirmed segmentation; 
    %   (2) confirmed entering cell; 
    %   (3) contours needs to evolve
    [Ps,newCellFrame,BMap,propagateIdx] = ConvertCellToContour(cellSemiGlobal,I,Options);

    % contour evolution
    if(numel(Ps)>0)
        [newPs, divisionIDX]=OpenActiveContour(I,Ps,BMap,Options);
    else
        newPs=[]; divisionIDX=[];
    end
    
    % update
    [cellFrame, cMat, maxID]=updateCellEachFrame(cellBlock,...
    newCellFrame, newPs, propagateIdx, matBlock{3}.Mat ,...
    [xdim,ydim], divisionIDX, maxID, Options);

    %DrawSegmentedArea2D(cellFrame{2},mat2gray(I),2);
    drawColorRegions(cellFrame{2}, [xdim,ydim], frameIdx ,cMap);
    
    saveas(gcf,[fpath,'sq',num2str(sq),'\track\img0',num2str(frameIdx+100),'.png'],'png');
    
    % write the first frame in the block to disk
    cellFrameTracked=cellFrame{1};
    save([fpath,'sq',num2str(sq),'/track_data/seg0',num2str(100+frameIdx),'.mat'],'cellFrameTracked');
    clear cellFrameTracked

    % update the block
    cellBlock{1}=cellFrame{2}; matBlock{1}.Mat = cMat;
    cellBlock{2}=cellFrame{3}; matBlock{2} = matBlock{3};
    for i=3:numFrameAhead+1
        cellBlock{i}=cellBlock{i+1};
        matBlock{i}=matBlock{i+1};
    end
    
    I1=I2;
    I2=I3;
    
    clear cellFrame cMat newPs divisionIDX Ps I newCellFrame cellSemiGlobal propagateIdx idxConsider
end

clear I1 I2 I3 
