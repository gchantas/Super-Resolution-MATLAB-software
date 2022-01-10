function result=applyBlockDiag(A,xblock,decFac)

%This to work should all elementary blocks matrices be -----> diagonal

if decFac==1

    result=A.*xblock;


else



    result{1}=addBlockVec (  applyBlockDiag(A{1,1},xblock{1},decFac/2), applyBlockDiag(A{1,2},xblock{2},decFac/2) ,decFac/2 ) ;
    

    result{2}=addBlockVec (  applyBlockDiag(A{2,1},xblock{1},decFac/2), applyBlockDiag(A{2,1},xblock{2},decFac/2) ,decFac/2 )  ;
    
      


end;