function [P, divisionIDX]=OpenActiveContour(I,P,BMap0,Options)

% [O,J]=OpenActiveContour(I,P,Options)
%  
% inputs,
%   I : An Image of type double preferable ranged [0..1]
%   P : List with coordinates descriping the rough contour N x 2
%   Options : A struct with all snake options
%   
% outputs,
%   O : List with coordinates of the final contour M x 2
%   J : Binary image with the segmented region
%
% options (general),
%  Option.Verbose : If true show important images, default false
%  Options.nPoints : Number of contour points, default 20
%  Options.Gamma : Time step, default 1
%  Options.Iterations : Number of iterations, default 10
%
% options (internal force)
%  Options.Alpha : Membrame energy  (first order), default 0.2
%  Options.Beta : Thin plate energy (second order), default 0.0

% options (Snake)
%  Options.Delta : stretching force due to length prior, default 1
%  Options.Kappa : Weight of repelling force, default 0.25

% Function is written by D.Kroon University of Twente (July 2010)
% Modified by Jianxu Chen (University of Notre Dame) at Jan 2015

% Process inputs
defaultoptions=struct('Verbose',false,'nPoints',20,'Alpha',0.2,'Beta',0.0,'Delta',1,...
    'Gamma',1,'Kappa',0.25,'Iterations',10);

if(~exist('Options','var')), 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))), 
        warning('snake:unknownoption','unknown options found');
    end
end

% If color image convert to grayscale
if(size(I,3)==3)
    I=rgb2gray(I); 
else
    I=mat2gray(I);
end

% Make the interal force matrix (smooth the contour)
S=SnakeInternalForceMatrix2D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);


% Make an uniform sampled contour description
P=InterpolateContourPoints2D(P,Options.nPoints,size(I));

%%%%% code for inspection %%%%%%%%%
% tmp=zeros(size(I));
% for i=1:1:numel(P)
%     tmp(P{i}.region>0)=i;
% end
% figure, imshow(tmp+1, rand(200,3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P=cellInfoUpdate(P,I);
if(Options.Verbose)
    figure(2), imshow(I), hold on; myHandle=drawContours(P,0,[],0);
end

% prepare the barrier map
BMap = processBMap(BMap0,Options.repelThresh);

% % Transform the Image into an External Energy Image
% Eext = ExternalForceImage2D(I,0.04, 2, 0.01 ,6);
%  
% % Make the external force (flow) field.
% Fx=ImageDerivatives2D(Eext,20,'x');
% Fy=ImageDerivatives2D(Eext,20,'y');
% Fext(:,:,1)=-Fx*2*20^2;
% Fext(:,:,2)=-Fy*2*20^2;
% 
% Fext=GVFOptimizeImageForces2D(Fext, 0.1, 5, 1.0);

for i=1:Options.Iterations    
    P=SnakeRegionUpdate(I,S,P,BMap,BMap0,Options.Gamma,Options.Kappa,Options.Delta,Options.repelThresh);  
    P=InterpolateContourPoints2D(P,Options.nPoints,size(I));
    P = cellInfoUpdate(P,I);
    
    % Show current contour
    if(Options.Verbose)
        myHandle=drawContours(P,i/Options.Iterations,myHandle,i);
   % else
   %     disp(['iteration: ',num2str(i)]);
    end
    
    if(stopCheck(P))
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check cell division
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

divisionIDX = checkDivision(P,I);
if(~isempty(divisionIDX))
    disp('division found!');
    %keyboard
    numNew = numel(divisionIDX);
    Pnew = cell(1,numNew*2);
    for i=1:1:numNew
        cid=divisionIDX(i);
        pp = P{cid}.pts;
        tarLen = 0.45*P{cid}.targetLength;
        p1 = pp(1:1:floor(0.3*Options.nPoints),:);
        p2 = pp(ceil(0.7*Options.nPoints):1:end,:);
        
        %%% copy infor in last frame %%%
        intensityLast = P{cid}.LastFrameIntensity;
        ptsLast=P{cid}.LastFramePts;
        lenLast=round(0.45*size(ptsLast));
        
        Pnew{2*i-1}=struct('pts',p1,'thickness',P{cid}.thickness,'length',0,...
        'targetLength',tarLen,'strip1',[],'strip2',[],'region',[],'intensity',[],...
        'normvec',[],'LastFrameIntensity',intensityLast,'LastFramePts',ptsLast(1:lenLast,:));
        Pnew{2*i}=struct('pts',p2,'thickness',P{cid}.thickness,'length',0,...
        'targetLength',tarLen,'strip1',[],'strip2',[],'region',[],'intensity',[],...
        'normvec',[],'LastFrameIntensity',intensityLast,'LastFramePts',ptsLast(end-lenLast+1:end,:));
    end
    Pnew=InterpolateContourPoints2D(Pnew,Options.nPoints,size(I));
    Pnew = cellInfoUpdate(Pnew,I);

    % update BMap
    J=DrawSegmentedArea2D(P,I,1);
    BMap1 = J | BMap0;
    for i=1:1:numNew
        BMap1(P{divisionIDX(i)}.region>0)=0;
    end
    BMap1 = bwareaopen(BMap1,50);
    
    BMap = processBMap(BMap1,Options.repelThresh);
    
    if(Options.Verbose)
         figure(2), imshow(I), hold on; myHandle=drawContours(Pnew,0,myHandle,0);
    end
    
    for i=1:Options.Iterations
        Pnew=SnakeRegionUpdate(I,S,Pnew,BMap,BMap1,Options.Gamma,Options.Kappa,Options.Delta,Options.repelThresh);
        Pnew=InterpolateContourPoints2D(Pnew,Options.nPoints,size(I));
        Pnew = cellInfoUpdate(Pnew,I);
        if(Options.Verbose)
            myHandle=drawContours(Pnew,i/Options.Iterations,myHandle,i);
        end
    end
    
    %%% merge divided contours back into P
    for i=1:1:numNew
        cid=divisionIDX(i);
        P{cid}=Pnew{2*i-1};
        P=cat(2,P,Pnew{2*i});
    end
end

    
    
%if(nargout>1)
%     J=DrawSegmentedArea2D(P,I,1);
%end

end

function flag=stopCheck(P)
    flag=true;
    for i=1:1:numel(P)
        %disp([abs(P{i}.length - P{i}.targetLength),0.05*P{i}.targetLength]);
        if(abs(P{i}.length - P{i}.targetLength)> 0.05*P{i}.targetLength)
            flag=false;
            break;
        end
    end
end

