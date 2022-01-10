# Super-Resolution-MATLAB-software
This is the code used to produce the results of the "Heavy tailed self-similarity modeling for Single Image Super Resolution". IEEE Transactions on Image Processing, 2020

G. Chantas, S. N. Nikolopoulos and I. Kompatsiaris, "Heavy-Tailed Self-Similarity Modeling for Single Image Super Resolution," in _IEEE Transactions on Image Processing_, vol. 30, pp. 838-852, 2021, doi: 10.1109/TIP.2020.3038521.


PDF: https://www.researchgate.net/publication/346922984_Heavy-Tailed_Self-Similarity_Modeling_for_Single_Image_Super_Resolution

See the usage examples of the RunExperiments.m file (used to produce the results in the paper) to see how the .m files are used:

pathHR: the path+name of the output file of the resulting high-res image
pathLR: the path+name of the input file of the low-res image
decFactor: the scaling factor which can be either 2 or 4. For decFactor=3 cases,there is a special function for this case (see below), so it is not passed as a parameter.


1. NonLocalPatchesFINDZ_YCbCr1_9_2019(pathHR,pathLR,imNumber,decFactor);
VBPS x2 or x4: This is the function implementing the standard Variational Bayes Super Resolution. decFactor is actually the scaling factor, which, here, can be **2 or 4**. 

2. NonLocalPatchesFINDZ_YCbCr_X3_22_7_2019(pathHR,pathLR,imNumber); 
VBPS x3: This is the function implementing the standard Variational Bayes Super Resolution for decFactor=3. 

3. SRviaDenoisingYCbCr22_1_2020(pathHR,pathSR,pathLR,decFactor)
t-VBPNLM x2 or x4: The function implementing the t-VBPNLM algorithm

4. SRviaDenoisingYCbCrGauss24_1_2020(pathHR,pathLR,decFactor);
VBPNLM x2 or x4: The function implementing the VBPNLM algorithm

5. SRviaDenoisingYCbCx3_r20_1_2020(pathHR,pathLR);
t-VBPNLM x3: The function implementing the t-VBPNLM (decFactor=3 in the function).

6. SRviaDenoisingYCbCx3Gauss_26_1_2020(pathHR,pathLR);
VBPNLM x3: The function implementing the t-VBPNLM (decFactor=3 in the function).
