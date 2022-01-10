%function BlockDiagMatCells = makeBlockDiagMatCells(Lhcong,Lh,Lq,Nx,Ny,decFac)
function  BlockDiagMatCells=makeDiagMatCellsOdd(Nx,Ny,decFac,iter)


if decFac==1
    %BlockDiagMatCells=Lhcong.*Lh+Lq;
    BlockDiagMatCells= zeros(Nx,Ny);%%ones(Nx,Ny);
else



    if mod(iter,2)==0
            Nx=Nx/3;
            Ny=Ny/3;
    end

decFac=decFac/3;
   
    BlockDiagMatCells{1,1}  =   makeDiagMatCellsOdd(   Nx,  Ny,  decFac, iter+1);
   
    BlockDiagMatCells{1,2}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
    
    BlockDiagMatCells{1,3}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
   
    BlockDiagMatCells{2,1}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
   
    BlockDiagMatCells{2,2}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
   
    BlockDiagMatCells{2,3}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
    
    BlockDiagMatCells{3,1}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
   
    BlockDiagMatCells{3,2}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
   
    BlockDiagMatCells{3,3}  =   makeDiagMatCellsOdd(   Nx,   Ny, decFac, iter+1);
   
   

   %BlockDiagMatCells{1,1}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)'  ,   Lh(1:Nx,1:Ny/2)', Lq(1:Nx,1:Ny/2)'  ,  Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{1,2}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)  ,   Lh(Nx/2+1:Nx,1:Ny), zeros(Nx/2,Ny/2), Nx/2,  Ny/2,   decFac/2);
   %BlockDiagMatCells{2,1}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1,1:Ny)  ,   Lh(1:Nx,1:Ny/2)   ,  zeros(Nx/2,Ny/2),   Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{2,2}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1:Nx,1:Ny)  ,  Lh(Nx/2+1:Nx,1:Ny),  Lq(Nx/2+1:Nx,1:Ny), Nx/2,  Ny/2,   decFac/2);

   %if makeAdd==1 && decFactor==2
   %makeAdd=0;
   %else
   %end

end