function diagonal=minusBlockDiagMat(A,decFac)

%This to work should all elementary blocks matrices be -----> diagonal

if decFac==1
    diagonal=-A;
else

   
    diagonal{1,1}=minusBlockDiagMat( A{1,1}  ,decFac/2)  ;

    diagonal{1,2}=minusBlockDiagMat( A{1,2}  ,decFac/2)  ;
    
    diagonal{2,1}=minusBlockDiagMat( A{2,1}  ,decFac/2)  ;
    
    diagonal{2,2}=minusBlockDiagMat( A{2,2}  ,decFac/2)  ;

    
end;


