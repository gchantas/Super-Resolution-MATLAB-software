function C=multiplyBlockDiagMat(A,B,decFac)

if decFac==1


    C=A.*B;


else
    
    decFac=decFac/3;
    
    C{1,1}  =   addBlockMat( multiplyBlockDiagMat( A{1,1}, B{1,1},   decFac) ,   multiplyBlockDiagMat( A{1,2}, B{2,1} ,   decFac), multiplyBlockDiagMat( A{1,3}, B{3,1} ,   decFac),decFac);
    C{1,2}  =   addBlockMat( multiplyBlockDiagMat( A{1,1}, B{1,2},   decFac) ,   multiplyBlockDiagMat( A{1,2}, B{2,2}  ,   decFac), multiplyBlockDiagMat( A{1,3}, B{3,2} ,   decFac),  decFac);
    C{1,3}  =   addBlockMat( multiplyBlockDiagMat( A{1,1}, B{1,3},   decFac) ,   multiplyBlockDiagMat( A{1,2}, B{2,3}  ,   decFac),multiplyBlockDiagMat( A{1,3}, B{3,3} ,   decFac),decFac);
    
    
    C{2,1}  =   addBlockMat( multiplyBlockDiagMat( A{2,1}, B{1,1},   decFac) ,   multiplyBlockDiagMat( A{2,2}, B{2,1}  ,   decFac),   multiplyBlockDiagMat( A{2,3}, B{3,1}, decFac) , decFac);
    C{2,2}  =   addBlockMat( multiplyBlockDiagMat( A{2,1}, B{1,2},   decFac) ,   multiplyBlockDiagMat( A{2,2}, B{2,2}  ,   decFac),   multiplyBlockDiagMat( A{2,3}, B{3,2}, decFac) ,decFac);
    C{2,3}  =   addBlockMat( multiplyBlockDiagMat( A{2,1}, B{1,3},   decFac) ,   multiplyBlockDiagMat( A{2,2}, B{2,3}  ,   decFac),  multiplyBlockDiagMat( A{2,3}, B{3,3}, decFac), decFac);
    
    C{2,1}  =   addBlockMat( multiplyBlockDiagMat( A{3,1}, B{1,1},   decFac) ,   multiplyBlockDiagMat( A{3,2}, B{2,1}  ,   decFac),   multiplyBlockDiagMat( A{3,3}, B{3,1}, decFac), decFac);
    C{2,2}  =   addBlockMat( multiplyBlockDiagMat( A{3,1}, B{1,2},   decFac) ,   multiplyBlockDiagMat( A{3,2}, B{2,2}  ,   decFac),   multiplyBlockDiagMat( A{3,3}, B{3,3}, decFac),  decFac);
    C{2,3}  =   addBlockMat( multiplyBlockDiagMat( A{3,1}, B{1,3},   decFac) ,   multiplyBlockDiagMat( A{3,2}, B{2,3}  ,   decFac),   multiplyBlockDiagMat( A{3,3}, B{3,3}, decFac),  decFac);

end