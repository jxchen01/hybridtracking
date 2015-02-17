function Options=setParameters(sq)

Options=struct();

Options.Verbose=false;

% Framework Parameters
Options.numFrameAhead = 3;

% Local EMD Matching Parameters
Options.candiRadius=20;
Options.bodyRatio=0.75;
Options.BoundThresh=5;
    
% Global EMD Matching Parameters
Options.searchDepth = min([3,Options.numFrameAhead+1]); % jump 3 frames in maximum
cK = zeros(1,Options.searchDepth-1);
cK(1) = Options.candiRadius;
for i=2:1:Options.searchDepth-1
    cK(i)= cK(i-1) + Options.candiRadius/2;
end
Options.candiRadiusK = cK;

% Contour Evolution Parameters
Options.Iterations=40;
Options.nPoints=20;
Options.ShrinkPixelNum = 12;

Options.Alpha=0.4;
Options.Beta=0.2;

if(sq==8)
    Options.minDij=5;
else
    Options.minDij=3;
end    

if(sq<10)
    Options.pixelCanSkip=7;
    Options.lengthCanSkip=10; % half contour length
    Options.repelThresh=6;
elseif(sq>10)
    Options.pixelCanSkip=10;
    Options.lengthCanSkip=14; % half contour length
    Options.repelThresh=8;
end

Options.leavingLength=20;

% Prior Information Parameters
Options.minBranch = 7; %%% used for pruning centerlines
