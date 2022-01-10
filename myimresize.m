function y=myimresize(x,Nx,Ny,filtDFTcoef,decFac,extras)



Nx=Nx+2*extras;
Ny=Ny+2*extras;


% hh1=zeros(Nx,Ny);
%  sigall=0.5;
%  
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




x_=padarray(x,[extras extras],'symmetric');



xfilt = real(ifft2( fft2(x_) .* (filtDFTcoef) ) );



y=xfilt(1+extras:decFac:Nx-extras,1+extras:decFac:Ny-extras);




