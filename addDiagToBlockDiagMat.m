function blockDiagMatCells = addDiagToBlockDiagMat(blockDiagMatCells,Lq,Nx,Ny,decFac)

if decFac==1
    size(Lh1)
    size(Lh2)
    size(Lq)
    blockDiagMatCells=Lh+Lq;
else


    blockDiagMatCells{1,1}  =   addDiagToBlockDiagMat( blockDiagMatCells{1,1} ,  Lq(1:Nx/2,1:Ny), Nx/2,  Ny/2,  decFac/2);




    blockDiagMatCells{2,2}  =   addDiagToBlockDiagMat( blockDiagMatCells{2,2},  Lq(Nx/2+1,N,1:Ny),  Nx/2,  Ny/2,  decFac/2);

end
    

