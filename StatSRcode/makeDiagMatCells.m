%function BlockDiagMatCells = makeBlockDiagMatCells(Lhcong,Lh,Lq,Nx,Ny,decFac)
function  BlockDiagMatCells=makeDiagMatCells(Diagonals,Nx,Ny,decFac,iter)


if decFac==1
    %BlockDiagMatCells=Lhcong.*Lh+Lq;
    BlockDiagMatCells=Diagonals;%%ones(Nx,Ny);
else

    if mod(iter,2)==0
        D1{1}= Diagonals(1:Nx/3,1:Ny/3); 
        D1{2}= Diagonals(Nx/3+1:2*Nx/3,1:Ny/3); 
        D1{3}= Diagonals(2*Nx/3+1:Nx,1:Ny/3); 
        
        D2{1}= Diagonals(1:Nx/3,Ny/3+1:2*Ny/3); 
        D2{2}= Diagonals(Nx/3+1:2*Nx/3,Ny/3+1:2*Ny/3); 
        D2{3}= Diagonals(2*Nx/3+1:Nx,Ny/3+1:2*Ny/3);

        D3{1}= Diagonals(1:Nx/3,2*Ny/3+1:Ny); 
        D3{2}= Diagonals(Nx/3+1:2*Nx/3,2*Ny/3+1:Ny);
        D3{3}= Diagonals(2*Nx/3+1:Nx,2*Ny/3+1:Ny); 
        Nx=Nx/3;
        Ny=Ny/3;
    else

        D1=Diagonals{1};
        D2=Diagonals{2};
        D3=Diagonals{3};
    end

   decFac=decFac/3;

   BlockDiagMatCells{1,1}  =   makeDiagMatCells(  D1, Nx,  Ny,  decFac,iter+1);
   BlockDiagMatCells{1,2}  =   makeDiagMatCellsOdd(     Nx,   Ny,   decFac,  iter+1);
   BlockDiagMatCells{1,3}  =   makeDiagMatCellsOdd(     Nx,   Ny,   decFac,  iter+1);

   BlockDiagMatCells{2,1}  =   makeDiagMatCellsOdd(     Nx,   Ny,   decFac,  iter+1);
   BlockDiagMatCells{2,2}  =   makeDiagMatCells(  D2,  Nx,   Ny,   decFac,iter+1);
   BlockDiagMatCells{2,3}  =   makeDiagMatCellsOdd(     Nx,   Ny,   decFac,  iter+1);

   BlockDiagMatCells{3,1}  =   makeDiagMatCellsOdd(     Nx,   Ny,   decFac,  iter+1);
   BlockDiagMatCells{3,2}  =   makeDiagMatCells( D3  ,   Nx,   Ny,   decFac,  iter+1);
   BlockDiagMatCells{3,3}  =   makeDiagMatCellsOdd(     Nx,   Ny,   decFac,  iter+1);

   %BlockDiagMatCells{1,1}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)'  ,   Lh(1:Nx,1:Ny/2)', Lq(1:Nx,1:Ny/2)'  ,  Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{1,2}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)  ,   Lh(Nx/2+1:Nx,1:Ny), zeros(Nx/2,Ny/2), Nx/2,  Ny/2,   decFac/2);
   %BlockDiagMatCells{2,1}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1,1:Ny)  ,   Lh(1:Nx,1:Ny/2)   ,  zeros(Nx/2,Ny/2),   Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{2,2}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1:Nx,1:Ny)  ,  Lh(Nx/2+1:Nx,1:Ny),  Lq(Nx/2+1:Nx,1:Ny), Nx/2,  Ny/2,   decFac/2);

   %if makeAdd==1 && decFactor==2
   %makeAdd=0;
   %else
   %end

end

