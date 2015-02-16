function drawColorRegions(cellList, sz, frameIdx, cMap)

cimg = uint16(zeros(sz));

for  i=1:1:numel(cellList)
   cimg(cellList{i}.seg>0) = cellList{i}.id;
end

figure(5); imshow(cimg, cMap), title(['frame: ',num2str(frameIdx)]), drawnow;