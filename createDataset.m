clear all
for imNumber=1:5
name = sprintf('/home/gchantas/Downloads/Set5/image_%03d.png',imNumber);

I_=uint8(imread(name));

[Nx,Ny,color]=size(I_);

dfac=2;
l2x=dfac*floor(Nx/dfac);
l2y=dfac*floor(Ny/dfac);

if color>1
I__=rgb2ycbcr(I_(1:l2x,1:l2y,:));


I=I__(1:l2x,1:l2y,1);

else
I=I_(1:l2x,1:l2y);
end
%I=imresize(I,1.0/dfac);
name2 = sprintf('/home/gchantas/Downloads/Set5x1_2/image_%03d.png',imNumber);

imwrite(I,name2);

end