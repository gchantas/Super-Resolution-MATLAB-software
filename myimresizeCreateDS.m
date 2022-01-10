function y=myimresizeCreateDS(x,decFac)

extras=3*4;%expand by $extras pixels the boundary of the image by replicating the bounds

[Nx,Ny]=size(x);

Nx=Nx+2*extras;
Ny=Ny+2*extras;


hh1=zeros(Nx,Ny);


sigall=1/decFac+0.01;

%Bicubic separable filter

    sigblurx=sigall;
    sigblury=sigall;

    for i=round(Nx/2)-12:round(Nx/2)+13
        
        for j=round(Ny/2)-12:round(Ny/2)+13
            if decFac==4
                hh1(i,j)=sigblurx*cubic(-(i-Nx/2+1/4)*sigblurx-1/16)*sigblury*cubic(-(j-Ny/2+1/4)*sigblury-1/16);
                %For decFactor=2
                
            elseif decFac==2
                hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);
            end
            
        end
    end


hh1=hh1/sum(sum(hh1));
hh1=fftshift(hh1);
filtDFTcoef=(fft2(hh1));

x_=padarray(x,[extras extras],'symmetric');



xfilt = real(ifft2( fft2(x_) .* (filtDFTcoef) ) );



y=xfilt(1+extras:decFac:Nx-extras,1+extras:decFac:Ny-extras);




