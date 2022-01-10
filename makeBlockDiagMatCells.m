%function BlockDiagMatCells = makeBlockDiagMatCells(Lhcong,Lh,Lq,Nx,Ny,decFac)
function BlockDiagMatCells = makeBlockDiagMatCells(Nx,Ny,decFac,iter)

%Iter should start with 1. It selects whether Nx and Ny will be divided
%(this is complicated due to two dimensional signal)

if decFac==1

    %BlockDiagMatCells=Lhcong.*Lh+Lq;
    BlockDiagMatCells=ones(Nx,Ny);

else
 
    
    if mod(iter,2)==1
        Nx=Nx/2;
        Ny=Ny/2;
    end
    BlockDiagMatCells{1,1}  =   makeBlockDiagMatCells(Nx,  Ny,  decFac/2,iter+1);
    BlockDiagMatCells{1,2}  =   makeBlockDiagMatCells(Nx,  Ny,  decFac/2,iter+1);
    BlockDiagMatCells{2,1}  =   makeBlockDiagMatCells(Nx,  Ny,  decFac/2,iter+1);
    BlockDiagMatCells{2,2}  =   makeBlockDiagMatCells(Nx,  Ny,  decFac/2,iter+1);

   %BlockDiagMatCells{1,1}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)'  ,   Lh(1:Nx,1:Ny/2)', Lq(1:Nx,1:Ny/2)'  ,  Ny/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{1,2}  =   makeBlockDiagMatCells(  Lhcong(1:Nx,1:Ny/2)  ,   Lh(Nx/2+1:Nx,1:Ny), zeros(Nx/2,Ny/2), Nx/2,  Ny/2,   decFac/2);
   %BlockDiagMatCells{2,1}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1,1:Ny)  ,   Lh(1:Nx,1:Ny/2)   ,  zeros(Nx/2,Ny/2),   Nx/2,   Ny/2,   decFac/2);
   %BlockDiagMatCells{2,2}  =   makeBlockDiagMatCells(  Lhcong(Nx/2+1:Nx,1:Ny)  ,  Lh(Nx/2+1:Nx,1:Ny),  Lq(Nx/2+1:Nx,1:Ny), Nx/2,  Ny/2,   decFac/2);
    
   %if makeAdd==1 && decFactor==2
   %makeAdd=0;
   %else
   %end

end