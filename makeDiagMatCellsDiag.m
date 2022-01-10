%function BlockDiagMatCells = makeBlockDiagMatCells(Lhcong,Lh,Lq,Nx,Ny,decFac)
function  BlockDiagMatCells=makeDiagMatCellsDiag(Diagonals,Nx,Ny,decFac)


if decFac==1

    %BlockDiagMatCells=Lhcong.*Lh+Lq;
    BlockDiagMatCells=ones(Nx,Ny);

else
 
    BlockDiagMatCells{1,1}  =   makeDiagMatCellsDiag(Ny/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{1,2}  =   makeDiagMatCells(Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{2,1}  =   makeDiagMatCells(Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{2,2}  =   makeDiagMatCellsDiag(Nx/2,  Ny/2,  decFac/2);

   %BlockDiagMatCells{1,1}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)'  ,   Lh(1:Nx,1:Ny/2)', Lq(1:Nx,1:Ny/2)'  ,  Ny/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{1,2}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)  ,   Lh(Nx/2+1:Nx,1:Ny), zeros(Nx/2,Ny/2), Nx/2,  Ny/2,   decFac/2);
   %BlockDiagMatCells{2,1}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1,1:Ny)  ,   Lh(1:Nx,1:Ny/2)   ,  zeros(Nx/2,Ny/2),   Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{2,2}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1:Nx,1:Ny)  ,  Lh(Nx/2+1:Nx,1:Ny),  Lq(Nx/2+1:Nx,1:Ny), Nx/2,  Ny/2,   decFac/2);
    
   %if makeAdd==1 && decFactor==2
   %makeAdd=0;
   %else
   %end

end