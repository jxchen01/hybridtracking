function R=regionRefine(I)

%%% fill small holes %%%
bg=imcomplement(I);
bgNew=bwareaopen(bg,6,4);
R=imcomplement(bgNew);

%%% remove small spur/bump on boundary %%%
nbrKernel=[0,1,0;1,0,1;0,1,0];
resp=conv2(single(R),single(nbrKernel),'same');
resp(~R)=0;
R(resp<2)=0;

