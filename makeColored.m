clear all
for imNumber=1:14
%name = sprintf('/home/gchantas/Downloads/Set5/image_%03d.png',imNumber);
dfac=4;

%name = sprintf('../Urban100/img_%03d.png',imNumber);
name = sprintf('../Set14/image_%03d.png',imNumber);
%name_ = sprintf('../Urban100x4_SRoutput/Urban100x4_%3d',imNumber);
name_ = sprintf('../Set14x4_SRoutput/Set14x4_%03d.png',imNumber);
%name2 = sprintf('/home/gchantas/Downloads/Set5x1_2/image_%03d.png',imNumber);
%name2 = sprintf('../Urban100x4color/Urban100x4_%03d.png',imNumber);
name2 = sprintf('../Set14x4color/Set14x4_%03d.png',imNumber);

I=im2double(imread(name));
I_=im2double(imread(name_));

[Nx,Ny,color]=size(I);


l2x=dfac*floor(Nx/dfac);
l2y=dfac*floor(Ny/dfac);

if color>1
I=rgb2ycbcr(I(1:l2x,1:l2y,:));

else
imwrite(I,name2,'png');
continue;
end



I2=zeros(l2x,l2y,3);

I2(:,:,1)=I_(:,:,1);
I2(:,:,2)=imresize(imresize(I(:,:,2),1/dfac),dfac);
%I=imresize(I,1.0/dfac);
I2(:,:,3)=imresize(imresize(I(:,:,3),1/dfac),dfac);



imwrite(ycbcr2rgb(I2),name2,'png');

end