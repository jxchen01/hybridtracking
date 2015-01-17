function repelMap=processBMap(BMap,repelThresh)

%%% key parameter %%%
%repelThresh=10;
%%%%%%%%%%%%%%%%%%%%%

BMap = BMap>0;

repelMap = bwdist(BMap);
repelMap(repelMap>repelThresh)=0;
idx=find(repelMap>0);
repelMap(idx)=1./(1+exp(2.*(repelMap(idx)-3)));
repelMap(BMap)=100;