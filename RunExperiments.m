%Run experiments
decFactor=4;
listPSNR=[];

listISNR=[];

listSSIM=[];

for imNumber=1:100
    imNumber
    tstart = tic;
pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);
%pathHR = sprintf('/home/gchantas/Downloads/Set14x1_%d/image_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Downloads/Set14x%d/image_%03d.png',decFactor,imNumber);

%pathHR = sprintf('/home/gchantas/Downloads/Set5x1_%d/image_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Downloads/Set5x%d/image_%03d.png',decFactor,imNumber);

    [l1 l2 l3]=NonLocalPatchesFINDZ_YCbCr2_1_2019(pathHR,pathLR,imNumber);
    listPSNR=[listPSNR,l1]; listISNR=[listISNR,l2]; listSSIM=[listSSIM,l3];

    
tend(imNumber)= toc(tstart);
%pause(15);
end;
