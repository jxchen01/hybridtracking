function Options=setParameters(xdim,ydim)

Options=struct();

% Framework Parameters
Options.numFrameAhead = 3;

% Local EMD Matching Parameters
Options.candiRadius=10;
Options.bodyRatio=0.75;
Options.BoundThresh=4;
halfDim = ceil(0.5 * max([xdim,ydim]));
if(mod(halfDim,200)>100  || halfDim<=200)
    Options.halfROIs=200;
else
    Options.halfROIs=200 + ceil( mod(halfDim,200) / floor(halfDim/200) );
end
    
% Global EMD Matching Parameters
Options.searchDepth = Options.numFrameAhead+1; % Match frame k with k+1, k+2, ..., k+numFrameAhead.
cK = zeros(1,Options.numFrameAhead);
for i=1:1:Options.numFrameAhead-1
    cK(i)= Options.candiRadius + 5*i;
end
Options.candiRadiusK = cK;

% Contour Evolution Parameters
Options.Verbose=true;
Options.Iterations=30;
Options.nPoints=20;
Options.ShrinkPixelNum = 10;
Options.lengthCanSkip=12;

% Prior Information Parameters
Options.minBranch = 7; %%% used for pruning centerlines
