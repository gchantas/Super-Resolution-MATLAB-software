function result=applyBlockDiag(A,xblock,decFac)

%This to work should all elementary blocks matrices be -----> diagonal

if decFac==1

    result=A.*xblock;


else


    decFac=decFac/2;

    result{1}  =   addBlockVec (  applyBlockDiag(A{1,1},xblock{1},decFac), applyBlockDiag(A{1,2},xblock{2},decFac) ,decFac ) ;
    

    result{2}  =   addBlockVec (  applyBlockDiag(A{2,1},xblock{1},decFac), applyBlockDiag(A{2,2},xblock{2},decFac) ,decFac )  ;
    
      


end;