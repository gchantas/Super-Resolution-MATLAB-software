function C=addBlockMat(A,B,decFac)

if decFac==1
   
    C=A+B;
else
    
    decFac=decFac/3;
    
    C{1,1} = addBlockMat( A{1,1}, B{1,1}, decFac);
    C{1,2} = addBlockMat( A{1,2}, B{1,2}, decFac);
    C{1,3} = addBlockMat( A{1,3}, B{1,3}, decFac);
    
    
    C{2,1} = addBlockMat( A{2,1}, B{2,1}, decFac);
    C{2,2} = addBlockMat( A{2,2}, B{2,2}, decFac);
    C{2,3} = addBlockMat( A{2,3}, B{2,3}, decFac);
    
    
    C{3,1} = addBlockMat( A{3,1}, B{3,1}, decFac);
    C{3,2} = addBlockMat( A{3,2}, B{3,2}, decFac);
    C{3,3} = addBlockMat( A{3,3}, B{3,3}, decFac);
    
    
    
end