%Run experiments
clear all
decFactor=4
listPSNR=[];

listISNR=[];

listSSIM=[];


%namesHR=dir('/home/gchantas/Downloads/Set14x1_2');
%namesLR=dir('/home/gchantas/Downloads/Set14x2');


%namesLR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/LR/LRBI/LR20x2');
%namesHR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/HR/HR20_x2');

namesHR = dir(sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/',decFactor));
namesLR= dir(sprintf('/home/gchantas/Documents/SR/Urban100x%d/',decFactor));
    
    
[num,dummy]=size(namesLR);


for imNumber=22:num


    imNumber
    tstart = tic;
    %pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
    %pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);
    pathLR = sprintf('%s/%s',namesLR(imNumber).folder,namesLR(imNumber).name)
    [~,ilen]=size(namesHR(imNumber).name);
    %pathSR = sprintf('/home/gchantas/Downloads/July2020runs/NDVI1/PNLMx2/%s.png',namesHR(imNumber).name(1:ilen-4));
    pathSR = sprintf('/home/gchantas/Downloads/July2020runs/Urban100/PNLMx4/%s.png',namesHR(imNumber).name(1:ilen-4));


    pathHR = sprintf('%s/%s',namesHR(imNumber).folder,namesHR(imNumber).name);


    % pathHR = sprintf('/home/gchantas/Downloads/Set5x1_%d/image_%03d.png',decFactor,imNumber);
    %pathLR= sprintf('/home/gchantas/Downloads/Set5x%d/image_%03d.png',decFactor,imNumber);

    %pathHR = sprintf('/home/gchantas/Documents/PIRM2018/Test/OriginalYCbCr/%d.png',imNumber);
    %pathLR = sprintf('/home/gchantas/Documents/PIRM2018/Test/4x_downsampledYCbCr/%d.png',imNumber);

  % [l1 l2 l3]=NonLocalPatchesFINDZ_YCbCr_X3_20_7_2020(pathHR,pathSR,pathLR,imNumber);
  % NonLocalPatchesFINDZ_YCbCr1_9_2019(pathHR,pathSR,pathLR,imNumber,decFactor)
   %NonLocalPatchesFINDZ_YCbCrFast24_6_2020(pathHR,pathLR,pathSR,imNumber,decFactor);

   % z=im2double(imread(pathLR));
   % y = wsdsr(z, decFactor);
    
   % imwrite(y,pathSR);

   % imwrite(z,pathSR);
    
 %[l1 l2 l3]=SRviaDenoisingYCbCrPNLM_7_7_2020(pathHR,pathSR,pathLR,decFactor);
%   [l1 l2 l3]=SRviaDenoisingYCbCrGauss24_1_2020(pathHR,pathLR,decFactor);
   % [l1 l2 l3]=SRviaDenoisingYCbCx3_r20_1_2020(pathHR,pathSR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCx3Gauss_26_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCr22_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCrPNLM_12_2_2020(pathHR,pathSR,pathLR,decFactor);

   [l1 l2 l3]=SRviaDenoisingYCbCrx4_21_2_2020(pathHR,pathLR,decFactor);%/x4

   % listPSNR=[listPSNR, l1];
   % listISNR=[listISNR, l2];
   % listSSIM=[listSSIM, l3];

    %mean(listPSNR)
    tend(imNumber) = toc(tstart);

    %pause(15);

end;                   