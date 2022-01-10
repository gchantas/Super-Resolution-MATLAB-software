function C=addBlockVec(A,B,decFac)

if decFac==1
   
    C=A+B;
else
    C{1}=addBlockVec( A{1}, B{1},decFac/2);
    C{2}=addBlockVec( A{2}, B{2},decFac/2);
end