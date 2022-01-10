function y1 = Amat_fast(x, iw, rw )

Nx=iw{5};
Ny=iw{6};
maxX=iw{7};
maxY=iw{8};

P1=size(rw');
P=P1(1);

decFactor=iw{2};
H1=iw{3};
extras=iw{4};
coord=zeros(P,2);
for k=1:2:P
coord(k,:)=rw{k};

coord(k,1)=mod(coord(k,1)+Nx/2,Nx)-Nx/2;%mod(k1,Nx);
coord(k,2)=mod(coord(k,2)+Ny/2,Ny)-Ny/2;%mod(l1,Ny);
end

%maxX=8;%max(abs(mod(coord(:,1)+3*Nx,Nx)))
%maxY=8;%max(abs(mod(coord(:,2)+3*Ny,Ny)))
%maxX=0;
%maxY=0;
nx1=(1:Nx)+maxX;
ny1=(1:Ny)+maxY;



x=reshape(x,Nx,Ny);




y_=myimresize(x,Nx,Ny,H1,decFactor,extras);

y=myimresizeTranspose(y_,Nx,Ny,H1,decFactor,extras);

%y=imresize(x,1/decFactor);
%  y_=real(ifft2( fft2(x) .* conj(H1)));
%  y=y_(1:decFactor:Nx,1:decFactor:Ny);
% %y=imresize(y,decFactor);
%  g1=gpuArray(zeros(Nx,Ny));
%  g1(1:decFactor:Nx,1:decFactor:Ny)=y;
%  y = real(ifft2( fft2(g1) .* H1));%imresize(g_r,decFactor);


x=padarray(x,[maxX maxY],'circular');
 

for k=1:2:P

    %temp=rw{k+1}.*real(  ifft2( rw{k}.*fft2(x)  ));

    %temp=rw{k+1}.* (x-circshift(x,rw{k}));
    temp=rw{k+1}(nx1,ny1).* (x(nx1,ny1)-x(nx1-coord(k,1), ny1-coord(k,2)));
    temp=padarray(temp,[maxX maxY],'circular');

    y = y  +   (temp(nx1,ny1)-temp(nx1+coord(k,1), ny1+coord(k,2)));

%    y = y+(temp - circshift(temp,-rw{k}) );
    
    
    %y=y+real(ifft2( conj(rw{k}).*fft2(temp) ));

end
%y=y(nx1,ny1);
%size(y)
y1=y(:);