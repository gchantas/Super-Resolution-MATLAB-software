function [x,alpha,ssigma,Inv_,C,xPseudo,InvC]=StatSR(x,g,alpha,ssigma,Nx,Ny,decFactor,f)

%%%%%%%%%%%%%%%%%%%%
%Generate Bicubic interpolator
hh1=zeros(Nx,Ny);
sigall=1/2;
 dfac=decFactor;
decFactor=decFactor^2/2;

%Bicubic separable filter

    sigblurx=sigall;
    sigblury=sigall;

    for i=round(Nx/2)-12:round(Nx/2)+13
        for j=round(Ny/2)-12:round(Ny/2)+13
      %      hh1(i,j)=sigblurx*cubic(-(i-Nx/2+1/4)*sigblurx-1/16)*sigblury*cubic(-(j-Ny/2+1/4)*sigblury-1/16);
                       
           %For decFactor=2
          hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);

        end
    end

hh1=hh1/sum(sum(hh1));
hh1=fftshift(hh1);
Hcubic=(fft2(hh1));

HS=abs(Hcubic).^2;
size(f)
size(HS)
Nx
Ny
decFactor
norm(myimresize(f,Nx,Ny,Hcubic,dfac,0)-imresize(f,1/dfac),'fro')^2/(Nx*Ny)


%%%%%%%%%%%%%%%%%%%%
%Generate Q-regularization operator

q=zeros(Nx,Ny);
q(1,1)=-4;
q(1,2)=1;
q(2,1)=1;
q(Nx,1)=1;
q(1,Ny)=1;
Q=fft2(q);

Qp=abs(Q).^2;

I=makeBlockDiagMatCells(Nx,Ny,decFactor*2,0);%%Block diagonal identity matrix
HcubBlock=makeDiagMatCells(Hcubic,Nx,Ny,decFactor*2,0);
HcubBlockConj=makeDiagMatCells(conj(Hcubic),Nx,Ny,decFactor*2,0);

C=multiplyBlockDiagMat(HcubBlockConj,multiplyBlockDiagMat(I,HcubBlock,decFactor*2),decFactor*2);


g1 = 4*myimresizeTranspose(g,Nx,Ny,Hcubic,dfac,0) ;


Z_=fft2(g1);

for iter=1:55


    Ia=makeDiagMatCells(  alpha*ssigma*Qp,Nx,Ny,decFactor*2,0);%%diagonal identity matrix (no blocks)
    C2=addBlockMat(C,Ia,decFactor*2);
    Inv=invertBlockMatrix(C2,decFactor*2);


    bVec=makeBlockVec(  Z_,Nx,Ny,decFactor,0);


    xblock=applyBlockDiag(Inv,bVec,decFactor*2);

    vec=deblockVec(xblock,Nx,Ny,decFactor*2,0);



    diagHHa  =   tageDiagBlockMat(Inv,Nx,Ny,decFactor*2,0);
    
   % HC=multiplyBlockDiagMat(HcubBlockConj,multiplyBlockDiagMat(Inv,HcubBlock,decFactor*2),decFactor*2);

   
    
    alpha= (Nx*Ny-1)  /   (  sum(sum(  real(diagHHa.*Qp   )))*ssigma    +   sum(sum( real(ifft2( Q.*vec ).^2 )))   );


   % Hvec=real(ifft2(vec.*(Hcubic)));

  %ssigma=(norm( Hvec(1:2:Nx,1:2:Ny)-g,'fro')^2  +   ssigma*sum(sum(real(  tageDiagBlockMat(multiplyBlockDiagMat(Inv,C,decFactor*2),Nx,Ny,decFactor*2,0))))  )/(Nx*Ny/decFactor^2)
end


alpha=alpha/2;

Ia=makeDiagMatCells(  alpha*ssigma*Qp,Nx,Ny,decFactor*2,0);%%diagonal identity matrix (no blocks)
C2=addBlockMat(C,Ia,decFactor*2);
Inv_=invertBlockMatrix(C2,decFactor*2);


bVec=makeBlockVec(  Z_,Nx,Ny,decFactor,0);


xblock=applyBlockDiag(Inv_,bVec,decFactor*2);

vec=deblockVec(xblock,Nx,Ny,decFactor*2,0);

x=real(ifft2(vec));

%     InvC=invertPseudoBlockMatrix(C,decFactor*2,10^(-8));
% 
%     xblock=applyBlockDiag(InvC,bVec,decFactor*2);
%     vec=deblockVec(xblock,Nx,Ny,decFactor*2,0);
% xPseudo=real(ifft2(vec));

