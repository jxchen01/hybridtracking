function flag=confirmEntry(cellEachFrame, idx)

frameIdxBase=1;
flag=false;

if(cellEachFrame{1}{idx}.length<5)
    return
end

numFrame = numel(cellEachFrame);

cid = cellEachFrame{1}{idx}.child;
if(isempty(cid))
    return
end

totalLen=0;

for kk=1:1:numel(cid)
    
    frameID = floor(cid(kk));
    cellID = uint16((cid(kk) - frameID)*1000);
    if(cellID==0)
        frameID = 2;
        cellID = uint16(cid(kk));
    else
        frameID = frameID-frameIdxBase;
    end
    
    if(frameID>numFrame)
        return;
    end
    
    totalLen = totalLen + cellEachFrame{frameID}{cellID}.length;
end

if(abs(cellEachFrame{1}{idx}.length - totalLen)<max([5,0.15*totalLen]) )
    flag=true;
end


