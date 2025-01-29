function newMask=removeDeadCellMask(mask,deadCellMask,timeRange)
% remove deadcell pixels from the mask input (nx,ny,nz,nt)
% input: original mask (nx,ny,nz,nt) & deadCellMask(nx,ny,nt)
% output: newMask (nx,ny,nz,nt)

if length(size(mask)) == 4
    nx=size(mask,1);
    ny=size(mask,2);
    nz=size(mask,3);
    nt=size(mask,4);
    newMask=false(nx,ny,nz,nt);
    for zz=1:nz
        for tt=timeRange
            mask1=mask(:,:,zz,tt);
            mask2=deadCellMask(:,:,tt);
            
            mask3=mask1-mask2;
            mask3(mask3==-1)=0;
            mask3=logical(mask3);
            newMask(:,:,zz,tt)=mask3;
        end
    end
elseif length(size(mask)) == 3
    nx=size(mask,1);
    ny=size(mask,2);
    nt=size(mask,3);
    newMask=false(nx,ny,nt);
    for tt=timeRange
        mask1=mask(:,:,tt);
        mask2=deadCellMask(:,:,tt);
        
        mask3=mask1-mask2;
        mask3(mask3==-1)=0;
        mask3=logical(mask3);
        newMask(:,:,tt)=mask3;
    end
end