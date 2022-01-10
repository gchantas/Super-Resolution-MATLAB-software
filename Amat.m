function y1 = Amat(x, iw, rw )

P1=size(rw');
P=P1(1);
N1=size(rw{2});
Nx=N1(1);
Ny=N1(2);
decFactor=iw{2};
H1=iw{3};
extras=iw{4};

x=reshape(x,Nx,Ny);
y_=myimresize(x,Nx,Ny,H1,decFactor,extras);

y=myimresizeTranspose(y_,Nx,Ny,H1,decFactor,extras);

%y=imresize(x,1/decFactor);
%  y_=real(ifft2( fft2(x) .* conj(H1)));
%  y=y_(1:decFactor:Nx,1:decFactor:Ny);
% %y=imresize(y,decFactor);
%  g1=gpuArray(zeros(Nx,Ny));
%  g1(1:decFactor:Nx,1:decFactor:Ny)=y;
%     
%  y = real(ifft2( fft2(g1) .* H1));%imresize(g_r,decFactor);
    
 


for k=1:2:P

    %temp=rw{k+1}.*real(  ifft2( rw{k}.*fft2(x)  ));

    temp=rw{k+1}.* (x-circshift(x,rw{k}));
   % y = y+(temp - imresize( imresize( circshift(temp,-rw{k}), 2), 0.5 ));

    y = y+(temp - circshift(temp,-rw{k}) );
    %y=y+real(ifft2( conj(rw{k}).*fft2(temp) ));
end

y1=y(:);