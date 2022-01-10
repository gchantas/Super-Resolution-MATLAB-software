# Super-Resolution-MATLAB-software
This is the code used to produce the results of the "Heavy tailed self-similarity modeling for Single Image Super Resolution". IEEE Transactions on Image Processing, 2020



See the usage examples of the RunExperiments.m file (used to produce the results in the paper) to see how the .m files are used:

        pathHR: the path+name of the output file of the resulting high-res image
        pathLR: the path+name of the input file of the low-res image
        decFactor: the scaling factor which can be either 2 or 4. For decFactor=3 cases,there is a special function for this case (see below), so it is not passed as a parameter.
        
        
    1. NonLocalPatchesFINDZ_YCbCr1_9_2019(pathHR,pathLR,imNumber,decFactor);
    This is the function implementing the standard Variational Bayes Super Resolution. decFactor is actually the scaling factor, which, here, can be **2 or 4**. 
    
    2. NonLocalPatchesFINDZ_YCbCr_X3_22_7_2019(pathHR,pathLR,imNumber); 
    This is the function implementing the standard Variational Bayes Super Resolution for decFactor=3. 
    
    3. SRviaDenoisingYCbCr14_1_2020(pathHR,pathLR,decFactor);
    The function implementing the 
    
    4. SRviaDenoisingYCbCrGauss24_1_2020(pathHR,pathLR,decFactor);
    5. SRviaDenoisingYCbCx3_r20_1_2020(pathHR,pathLR,decFactor);
    6. SRviaDenoisingYCbCx3Gauss_26_1_2020(pathHR,pathLR,decFactor);
    7. SRviaDenoisingYCbCr22_1_2020(pathHR,pathSR,pathLR,decFactor);
    8. SRviaDenoisingYCbCrPNLM_RandomD_11_6_2020(pathHR,pathLR,decFactor);
    9. SRviaDenoisingYCbCrx4_21_2_2020(pathHR,pathLR,decFactor);%/x4
