function C=addBlockMat(A,B,decFac)

if decFac==1
   
    C=A+B;
else
    C{1,1}=addBlockMat( A{1,1}, B{1,1}, decFac/2);
    C{1,2}=addBlockMat( A{1,2}, B{1,2}, decFac/2);
    C{2,1}=addBlockMat( A{2,1}, B{2,1}, decFac/2);
    C{2,2}=addBlockMat( A{2,2}, B{2,2}, decFac/2);
end