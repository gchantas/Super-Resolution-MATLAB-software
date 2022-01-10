function y=myimresizeTranspose(x,Nx,Ny,filtDFTcoef,decFac,extras)

Nx=Nx+2*extras;
Ny=Ny+2*extras;

x_=gpuArray2(zeros(Nx,Ny));


x_(1:decFac:Nx,1:decFac:Ny) = gpuArray2(  padarray(gather(x)  ,   [extras/decFac extras/decFac],'symmetric')   );



% hh1=zeros(Nx,Ny);
%  sigall=0.5;
%     
%     sigblurx=sigall;
%     sigblury=sigall;
%     for i=Nx/2-14:Nx/2+15
%         for j=Ny/2-14:Ny/2+15
%             hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);
%         end
%     end
%     
% hh1=hh1/sum(sum(hh1));
% hh1=fftshift(hh1);
% filtDFTcoef=(fft2(hh1));



%x_(1+extras:decFac:Nx-extras,1+extras:decFac:Ny-extras)=x;



%x_=padarray(x_,[extras extras],'replicate');




y_= real(ifft2( fft2(x_) .* conj(filtDFTcoef ) ));
y=y_(1+extras:Nx-extras,1+extras:Ny-extras);







