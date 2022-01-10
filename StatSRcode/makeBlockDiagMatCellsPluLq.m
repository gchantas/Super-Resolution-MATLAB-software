function BlockDiagMatCells = makeBlockDiagMatCellsPlusLq(Lh1,Lh2,Lq,Nx,Ny,decFac)

if decFac==1
    BlockDiagMatCells=Lh1.*Lh2+Lq;
else
    BlockDiagMatCells{1,1}  =   makeBlockDiagMatCellsPlusLq(  Lh1(1:Nx/2,1:Ny/2),  Lh2(1:Nx/2,1:Ny/2),  Lq, Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{1,2}  =   makeBlockDiagMatCells(  Lh1(1:Nx/2,Ny/2+1:Ny),  Lh2(1:Nx/2,Ny/2+1:Ny),  Lq, Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{2,1}  =   makeBlockDiagMatCells(  Lh1(Nx/2+1:Nx,Ny/2+1:Ny),  Lh2(Nx/2+1:Nx,1:Ny/2), Lq,  Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{2,2}  =   makeBlockDiagMatCellsPlusLq(  Lh1(Nx/2+1:Nx,Ny/2+1:Ny),  Lh2(Nx/2+1:Nx,Ny/2+1:Ny),  Lq, Nx/2,  Ny/2,  decFac/2);
    
    
end
    

