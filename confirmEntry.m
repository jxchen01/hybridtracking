function flag=confirmEntry(cellEachFrame, idx)

frameIdxBase=1;

numFrame = numel(cellEachFrame);


cid = cellEachFrame{1}{idx}.child;
frameID = floor(cid);
cellID = uint16((cid - frameID)*1000);
frameID = frameID-frameIdxBase;

flag=false;
if(frameID>numFrame)
    return;
end

if(abs(cellEachFrame{1}{idx}.length - cellEachFrame{frameID}{cellID}.length)...
        <max([5,0.15*cellEachFrame{frameID}{cellID}.length]) )
    flag=true;
end
    