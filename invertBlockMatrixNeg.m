function diagonalInverse=minusBlockDiagMat(A,decFac)

%This to work should all elementary blocks matrices be -----> diagonal

if decFac==1
    diagonalInverse=-1./A;
else

   
    diagonalInverse{1,1}=invertBlockMatrix(  subtractBlockDiagMat(A11, multiplyBlockDiagMat( A12, multiplyBlockDiagMat(A22inv, A21,decFac/2 ), decFac/2 ))  ,decFac/2)  ;
    
    
    temp1= invertBlockMatrix( subtractBlockDiagMat( multiplyBlockDiagMat(  multiplyBlockDiagMat( A21, A11inv, decFac/2) , A12 , decFac-1) , A22, decFac-1), decFac/2);
    diagonalInverse{1,2}=invertBlockMatrix(  multiplyBlockDiagMat(  multiplyBlockDiagMat( A11, A12, decFac/2), temp1 ,decFac/2), decFac/2) ;

    diagonalInverse{2,1}=diagonalInverse{2,1};%invertBlockMatrix(  subtractBlockDiagMat(B11, multiplyBlockDiagMat( B12, multiplyBlockDiagMat(A12inv, B21,decFac-1 ), decFac-1 ))  ,decFac-1)  ;

    diagonalInverse{2,2}=invertBlockMatrixNeg( A12 , subtractBlockDiagMat( multiplyBlockDiagMat(  multiplyBlockDiagMat( A21, A11inv, decFac/2) ,  decFac/2) , A22, decFac-1), decFac/2);

    
end;


