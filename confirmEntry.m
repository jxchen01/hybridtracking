function flag=confirmEntry(cellEachFrame, idx, LengthSkip,checkType)

%%%% check type:
% 1. entering from boundary
% 2. re-appearing, not close to boundary
frameIdxBase=1;
flag=false;

if(cellEachFrame{1}{idx}.length<LengthSkip)
    return
end

numFrame = numel(cellEachFrame);

cid = cellEachFrame{1}{idx}.child;
if(isempty(cid))
    return
end

if(checkType==1)

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
        return;
    end

elseif(checkType==2)
    
    if(cellEachFrame{1}{idx}.length<LengthSkip*2.5)
        return
    end

    if(numel(cid)~=1)
        return;
    end
    
    frameID = floor(cid);
    cellID = uint16((cid - frameID)*1000);
    if(cellID==0)
        frameID = 2;
        cellID = uint16(cid);
    else
        frameID = frameID-frameIdxBase;
    end
    
    if(abs(cellEachFrame{frameID}{cellID}.length - cellEachFrame{1}{idx}.length)...
            > max([5,0.15*cellEachFrame{frameID}{cellID}.length]))
        return;
    end
    
    ccid = cellEachFrame{frameID}{cellID}.child;
    if(isempty(ccid))
        return
    end
    
    totalLen=0;
    
    for kk=1:1:numel(ccid)
        
        cframeID = floor(ccid(kk));
        ccellID = uint16((ccid(kk) - cframeID)*1000);
        if(ccellID==0)
            cframeID = frameID+1;
            ccellID = uint16(ccid(kk));
        else
            cframeID = cframeID-frameIdxBase;
        end
        
        if(cframeID>numFrame)
            return;
        end
        
        totalLen = totalLen + cellEachFrame{cframeID}{ccellID}.length;
    end
    
    if(abs(cellEachFrame{frameID}{cellID}.length - totalLen)<max([5,0.15*totalLen]) )
        flag=true;
        return;
    end
    
        
else
    error('error in check type in confirmEntry');
end
    
    


