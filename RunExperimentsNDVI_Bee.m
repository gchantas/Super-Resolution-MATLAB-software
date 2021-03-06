%Run experiments
clear all
decFactor=3
listPSNR=[];

listISNR=[];

listSSIM=[];

namesLR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/LR/LRBI/LR20x3');
namesHR=dir('/home/gchantas/Downloads/DRLN-master/TestCode/HR/HR20_x3');
[num,dummy]=size(namesLR);


for imNumber=3:num


    imNumber
    tstart = tic;
    %pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
    %pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);
    pathLR = sprintf('%s/%s',namesLR(imNumber).folder,namesLR(imNumber).name)
    [~,ilen]=size(namesHR(imNumber).name);
    pathSR = sprintf('/home/gchantas/Downloads/DRLN-master/TestCode/SR/BI/DRLN_Set5/Set5/VBPSfullx3/%s.png',namesHR(imNumber).name(1:ilen-4));


    pathHR = sprintf('%s/%s',namesHR(imNumber).folder,namesHR(imNumber).name);


    % pathHR = sprintf('/home/gchantas/Downloads/Set5x1_%d/image_%03d.png',decFactor,imNumber);
    %pathLR= sprintf('/home/gchantas/Downloads/Set5x%d/image_%03d.png',decFactor,imNumber);

    %pathHR = sprintf('/home/gchantas/Documents/PIRM2018/Test/OriginalYCbCr/%d.png',imNumber);
    %pathLR = sprintf('/home/gchantas/Documents/PIRM2018/Test/4x_downsampledYCbCr/%d.png',imNumber);

 [l1 l2 l3]=NonLocalPatchesFINDZ_YCbCr_X3_25_7_2020(pathHR,pathSR,pathLR,imNumber);
    %  NonLocalPatchesFINDZ_YCbCr1_9_2019(pathHR,pathSR,pathLR,imNumber,decFactor);

   % z=im2double(imread(pathLR));
   % y = wsdsr(z, decFactor);
    
   % imwrite(y,pathSR);

   % imwrite(z,pathSR);
    
 %[l1 l2 l3]=SRviaDenoisingYCbCrPNLM_7_7_2020(pathHR,pathSR,pathLR,decFactor);
%   [l1 l2 l3]=SRviaDenoisingYCbCrGauss24_1_2020(pathHR,pathLR,decFactor);
%    [l1 l2 l3]=SRviaDenoisingYCbCx3_r20_1_2020(pathHR,pathSR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCx3Gauss_26_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCr22_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCrPNLM_12_2_2020(pathHR,pathLR,decFactor);

    %[l1 l2 l3]=SRviaDenoisingYCbCrx4_21_2_2020(pathHR,pathLR,decFactor);%/x4

   % listPSNR=[listPSNR, l1];
   % listISNR=[listISNR, l2];
   % listSSIM=[listSSIM, l3];

    %mean(listPSNR)
    tend(imNumber) = toc(tstart);

    %pause(15);

end;                   