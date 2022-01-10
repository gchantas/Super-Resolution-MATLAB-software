%function BlockDiagMatCells = makeBlockDiagMatCells(Lhcong,Lh,Lq,Nx,Ny,decFac)
function  BlockDiagMatCells=makeDiagMatCellsOdd2(Diagonals,Nx,Ny,decFac)


if decFac==1

    %BlockDiagMatCells=Lhcong.*Lh+Lq;
    BlockDiagMatCells=Diagonals;%%ones(Nx,Ny);

else
    
    
    
    BlockDiagMatCells{1,1}  =   makeDiagMatCellsOdd(Diagonals(1:Nx,1:Ny/2), Nx/1,  Ny/2,  decFac/2);
    BlockDiagMatCells{1,2}  =   makeDiagMatCells(zeros(Nx/2,Ny/2), Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{2,1}  =   makeDiagMatCells(zeros(Nx/2,Ny/2), Nx/2,  Ny/2,  decFac/2);
    BlockDiagMatCells{2,2}  =   makeDiagMatCellsOdd(Diagonals(1:Nx,Ny:Ny/2+1), Nx/1,  Ny/2,  decFac/2);

   %BlockDiagMatCells{1,1}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)'  ,   Lh(1:Nx,1:Ny/2)', Lq(1:Nx,1:Ny/2)'  ,  Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{1,2}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)  ,   Lh(Nx/2+1:Nx,1:Ny), zeros(Nx/2,Ny/2), Nx/2,  Ny/2,   decFac/2);
   %BlockDiagMatCells{2,1}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1,1:Ny)  ,   Lh(1:Nx,1:Ny/2)   ,  zeros(Nx/2,Ny/2),   Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{2,2}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1:Nx,1:Ny)  ,  Lh(Nx/2+1:Nx,1:Ny),  Lq(Nx/2+1:Nx,1:Ny), Nx/2,  Ny/2,   decFac/2);
    
   %if makeAdd==1 && decFactor==2
   %makeAdd=0;
   %else
   %end

end