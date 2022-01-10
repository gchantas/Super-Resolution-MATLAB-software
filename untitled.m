function y= expandim(x,Nx,Ny,bound)


x1=zeros(Nx+d,Ny+d);

Imask=zeros(size(x1));

Imask(1:Nx,1:Ny)=ones(Nx,Ny);

y=zeros(Nx,Ny);%returned image



x1(1:Nx,1:Ny)=x;


