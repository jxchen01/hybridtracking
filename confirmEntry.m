function flag=confirmEntry(cellEachFrame, idx)

frameIdxBase=1;
flag=false;

numFrame = numel(cellEachFrame);

cid = cellEachFrame{1}{idx}.child;
if(isempty(cid))
    return
end

frameID = floor(cid);
cellID = uint16((cid - frameID)*1000);
if(cellID==0)
    frameID = 2;
    cellID = uint16(cid);
else
    frameID = frameID-frameIdxBase;
end

if(frameID>numFrame)
    return;
end

if(abs(cellEachFrame{1}{idx}.length - cellEachFrame{frameID}{cellID}.length)...
        <max([5,0.15*cellEachFrame{frameID}{cellID}.length]) )
    flag=true;
end

