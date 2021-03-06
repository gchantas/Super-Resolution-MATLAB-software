function y1 = Amat_Q(x, iw, rw )

% Input:
%     iw{1}=1;
%     iw{2}=decFactor;
%     iw{3}=Hcubic;
%     iw{4}=extras;
%     iw{5}=Nx;
%     iw{6}=Ny;
%     iw{7}=maxX;
%     iw{8}=maxY;
% 
%     iw{9}=Inv;
%iw{10}=     ZeroIndX;

%iw{11}=ZeroIndY;
decFactor=iw{2};
H1=iw{3};
extras=iw{4};
Nx=iw{5};
Ny=iw{6};
x=reshape(x,Nx,Ny);
HDDH=iw{7};

P1=size(rw');
P=P1(1);

% 
X=fft2(x);
% 
% bVec=makeBlockVec( X,Nx,Ny,decFactor,0);
% 
% 
% xblock=applyBlockDiag(HDDH,bVec,decFactor*2);
% 
% Y=deblockVec(xblock,Nx,Ny,decFactor*2,0)/decFactor^2;

   % y=real(ifft2(vec))/decFactor^2;
y_=myimresize(x,Nx,Ny,H1,decFactor,extras);

y=myimresizeTranspose(y_,Nx,Ny,H1,decFactor,extras);
Y=fft2(y);
%y=imresize(x,1/decFactor);
%  y_=real(ifft2( fft2(x) .* conj(H1)));
%  y=y_(1:decFactor:Nx,1:decFactor:Ny);
% %y=imresize(y,decFactor);
%  g1=gpuArray(zeros(Nx,Ny));
%  g1(1:decFactor:Nx,1:decFactor:Ny)=y;
%  y = real(ifft2( fft2(g1) .* H1));%imresize(g_r,decFactor);

for k=1:2:P



    temp=rw{k+1}.*real(  ifft2( rw{k}.*X  ));

   Y=Y+conj(rw{k}).*fft2(temp) ;

end
y=real(ifft2(Y));
%y=y(nx1,ny1);
%size(y)
y1=y(:);