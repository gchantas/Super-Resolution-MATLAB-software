function x=RunSRcommands(pathLR,imNumber)

decFactor=4;


tstart = tic;
x=imcrop(im2double(imread(pathLR)));

imwrite(x,pathLR);
%pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);
%pathHR = sprintf('/home/gchantas/Downloads/Set14x1_%d/image_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Downloads/Set14x%d/image_%03d.png',decFactor,imNumber);
%pathHR = sprintf('/home/gchantas/Downloads/Set5x1_%d/image_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Downloads/Set5x%d/image_%03d.png',decFactor,imNumber);

pathHR = sprintf('/home/gchantas/Documents/AppDemo/hr_%d.png',imNumber);
pathLR = sprintf('/home/gchantas/Documents/AppDemo/original_%d.png',imNumber);

%[l1 l2 l3]=NonLocalPatchesFINDZ_YCbCr_X3_22_7_2019(pathHR,pathLR,imNumber);
[l1 l2 l3]=NonLocalPatchesFINDZ_YCbCr1_9_2019(pathHR,pathLR,imNumber,decFactor);
listPSNR=[listPSNR,l1]; listISNR=[listISNR,l2]; listSSIM=[listSSIM,l3];

mean(listPSNR)
tend= toc(tstart);


%pause(15);