function C=subtractBlockDiagMat(A,B,decFac)

if decFac==1
    C=A-B;
else
    decFac=decFac/2;
    C{1,1}=subtractBlockDiagMat( A{1,1}, B{1,1},decFac);
    C{1,2}=subtractBlockDiagMat( A{1,2}, B{1,2},decFac);
    C{2,1}=subtractBlockDiagMat( A{2,1}, B{2,1},decFac);
    C{2,2}=subtractBlockDiagMat( A{2,2}, B{2,2},decFac);
end