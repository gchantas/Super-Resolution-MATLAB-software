% clear all
% 
% for imNumber=1:100
% 
%     name = sprintf('/home/gchantas/Documents/PIRM2018/Test/4x_downsampled/%d.png',200+imNumber);
% 
%     I_=uint8(imread(name));
% 
%     [Nx,Ny,color]=size(I_);
% 
%     dfac=1;
%     l2x=dfac*floor(Nx/dfac);
%     l2y=dfac*floor(Ny/dfac);
% 
%     if color>1
%         I__=rgb2ycbcr(I_(1:l2x,1:l2y,:));
%         I=I__(1:l2x,1:l2y,1);
%     else
%         I=I_(1:l2x,1:l2y);
%     end
% 
%     %I=myimresizeCreateDS(I,dfac)/255;
%     name2 = sprintf('/home/gchantas/Documents/PIRM2018/Test/4x_downsampledYCbCr/%d.png',imNumber);
% 
%     imwrite(I,name2);
% 
% end


clear all
namesLR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/LR/LRBI/LR20x3');
namesHRx3=dir('/home/gchantas/Downloads/DRLN-master/TestCode/HR/HR20_x3');
namesHR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/HR/HR22');
%namesHR=dir('/home/gchantas/Documents/SR/Urban100x1_2');
%namesSR2=dir('/home/gchantas/Downloads/DRLN-master/TestCode/LR/LRBI/Set5/x4');
%namesSR2=dir('/home/gchantas/Downloads/DRLN-master/TestCode/SR/BI/DRLN_Urban100/Urban100/x2');
%namesSR=dir('/home/gchantas/Documents/SR/wsdsr-v1.0.2/SR20x3/');
%namesSR=dir('/home/gchantas/Documents/Urban100x2June2019');
%namesSR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/SR/BI/DRLN_Set5/SR22/');
%namesSR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/SR/BI/DRLN_Set5/SR22_fast');

[num,dummy]=size(namesLR);


for imN=1:num-2
    namesHR(imN+2).name

f_=im2double(imread( sprintf('%s/%s',namesHR(imN+2).folder, namesHR(imN+2).name )  ));

f(1:1023,1:1023)=f_(1:1023,1:1023);
imwrite(f,sprintf('%s/%s',namesHRx3(1).folder, namesLR(imN+2).name));

g2=imresize(f,1/3);
imwrite(g2,sprintf('%s/%s',namesLR(1).folder, namesLR(imN+2).name));
end