function BlockDiagMatCells = addDiagToBlockDiagMat(Lh,Lq,Nx,Ny,decFac)

if decFac==1
    size(Lh1)
    size(Lh2)
    size(Lq)
    BlockDiagMatCells=Lh+Lq;
else
    BlockDiagMatCells{1,1}  =   addDiagToBlockDiagMat(  Lh1(1:Nx/2,1:Ny/2),  Lh2(1:Nx/2,1:Ny/2),  Lq(1:Nx/2,1:Ny/2), Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{2,2}  =   addDiagToBlockDiagMat(  Lh1(Nx/2+1:Nx,Ny/2+1:Ny),  Lh2(Nx/2+1:Nx,Ny/2+1:Ny),  Lq(Nx/2+1:Nx,Ny/2+1:Ny), Nx/2,  Ny/2,  decFac/2);
    
    
end
    

