%Run experiments
decFactor=2
listPSNR=[];

listISNR=[];

listSSIM=[];

for imNumber=1:1


    imNumber
    tstart = tic;
    %pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
    %pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);
    pathHR = sprintf('D:\\ImagingData\\bee_HR_1.png');
    pathLR= sprintf('D:\\ImagingData\\bee_LR_1.png');

    %pathHR = sprintf('/home/gchantas/Downloads/Set5x1_%d/image_%03d.png',decFactor,imNumber);
    %pathLR= sprintf('/home/gchantas/Downloads/Set5x%d/image_%03d.png',decFactor,imNumber);

    %pathHR = sprintf('/home/gchantas/Documents/PIRM2018/Test/OriginalYCbCr/%d.png',imNumber);
    %pathLR = sprintf('/home/gchantas/Documents/PIRM2018/Test/4x_downsampledYCbCr/%d.png',imNumber);

    %[l1 l2 l3]=NonLocalPatchesFINDZ_YCbCr_X3_22_7_2019(pathHR,pathLR,imNumber);
    [l1 l2 l3]=NonLocalPatchesFINDZ_YCbCr18_6_2020(pathHR,pathLR,imNumber,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCr14_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCrGauss24_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCx3_r20_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCx3Gauss_26_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCr22_1_2020(pathHR,pathLR,decFactor);
    %[l1 l2 l3]=SRviaDenoisingYCbCrPNLM_12_2_2020(pathHR,pathLR,decFactor);

    %[l1 l2 l3]=SRviaDenoisingYCbCrx4_21_2_2020(pathHR,pathLR,decFactor);%/x4

    listPSNR=[listPSNR, l1];
    listISNR=[listISNR, l2];
    listSSIM=[listSSIM, l3];

    mean(listPSNR)
    tend(imNumber) = toc(tstart);

    %pause(15);

end;