function Options=setParameters(xdim,ydim)

Options=struct();

% Framework Parameters
Options.numFrameAhead = 3;

% Local EMD Matching Parameters
Options.candiRadius=20;
Options.bodyRatio=0.75;
Options.BoundThresh=4;
% halfDim = ceil(0.5 * max([xdim,ydim]));
% if(mod(halfDim,200)>100  || halfDim<=200)
%     Options.halfROIs=200;
% else
%     Options.halfROIs=200 + ceil( mod(halfDim,200) / floor(halfDim/200) );
% end
    
% Global EMD Matching Parameters
Options.searchDepth = min([3,Options.numFrameAhead+1]); % jump 3 frames in maximum
cK = zeros(1,Options.searchDepth-1);
cK(1) = Options.candiRadius;
for i=2:1:Options.searchDepth-1
    cK(i)= cK(i-1) + Options.candiRadius/2;
end
Options.candiRadiusK = cK;

% Contour Evolution Parameters
Options.Verbose=false;

Options.Iterations=30;
Options.nPoints=20;
Options.ShrinkPixelNum = 12;
Options.repelThresh=6;
Options.Alpha=0.4;
Options.Beta=0.2;

Options.lengthCanSkip=7;
Options.leavingLength=20;

Options.maxNormMove=3;
Options.maxTangMove=20;

% Prior Information Parameters
Options.minBranch = 7; %%% used for pruning centerlines
