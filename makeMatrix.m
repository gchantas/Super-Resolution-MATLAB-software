
clear all
Nx=1024;
Ny=1024;
df=2;

itnum=1000;

Hmat=zeros(Nx/2,Ny/2);

for iter1=1:itnum
    randn('state',iter1);
    %d=gpuArray(zeros(Nx,Ny));d(Nx/2+iter1,Ny/2+iter1)=1;
    d{iter1}=randn(Nx,Ny);
    e{iter1}=zeros(Nx,Ny);
e{iter1}=(imresize((d{iter1}),1/df));

%Hmat=Hmat+e.^2/itnum;

end

imNumber=801;

name = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_HR/0%d.png',imNumber);


f_= im2double(rgb2gray(imread(name))) ;

[Nx, Ny]=size(f_);
 Nx=min(1024,Nx);
 Ny=min(1024,Ny);
%f=f_(100:Nx+99,100:Nx+99);
f=(f_(1:Nx,1:Ny));


f2=zeros(Nx/df,Ny/df);


for iter1=1:itnum
    %randn('state',iter1);
    %d=gpuArray(zeros(Nx,Ny));
   % d(Nx/2+iter1,Ny/2+iter1)=1;
     %   d=gpuArray(sign(randn(Nx,Ny)));

    f2=sum(sum(d{iter1}.*f))*e{iter1}/itnum;

%Hmat=Hmat+e.^2/itnum;

end