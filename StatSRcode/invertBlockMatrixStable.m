function diagonalInverse=invertBlockMatrixStable(A,decFac)

%This to work should all elementary blocks matrices be -----> diagonal

if decFac==1

    diagonalInverse=1./A;


else

    A11=A{1,1};
    A22=A{2,2};
    A12=A{1,2};
    A21=A{2,1};

    A11inv=invertBlockMatrix(A11,decFac/2);
    %A12inv=invertBlockMatrix(B12,decFac-1);
    %A21inv=A12inv;%invertBlockMatrix(B21,decFac-1);
    A22inv=invertBlockMatrix(A22,decFac/2);



%     diagonalInverse{1,1}=invertBlockMatrix(  subtractBlockDiagMat(A11, multiplyBlockDiagMat( A12, multiplyBlockDiagMat(A22inv, A21, decFac/2 ), decFac/2 ) ,decFac/2 )  ,decFac/2);
%     temp1= invertBlockMatrix( subtractBlockDiagMat( multiplyBlockDiagMat(  multiplyBlockDiagMat( A21, A11inv, decFac/2) , A12 , decFac/2) , A22, decFac/2), decFac/2);
%     diagonalInverse{1,2}=invertBlockMatrix(  multiplyBlockDiagMat(  multiplyBlockDiagMat( A11inv, A12, decFac/2), temp1 ,decFac/2), decFac/2);
%     diagonalInverse{2,1}=diagonalInverse{1,2};%invertBlockMatrix(  subtractBlockDiagMat(B11, multiplyBlockDiagMat( B12, multiplyBlockDiagMat(A12inv, B21,decFac-1 ), decFac-1 ))  ,decFac-1)  ;
%     diagonalInverse{2,2}=minusBlockDiagMat(temp1,decFac/2);

  
    temp1= invertBlockMatrix( minus( multiplyBlockDiagMat( A11inv,  decFac/2) , A12 , decFac/2) , A22, decFac/2), decFac/2);


    diagonalInverse{1,1}=multiplyBlockDiagMat( A11inv, A12, decFac/2);
    
    
   
    
    diagonalInverse{1,2}=invertBlockMatrix(  multiplyBlockDiagMat(  multiplyBlockDiagMat( A11inv, A12, decFac/2), temp1 ,decFac/2), decFac/2);
    diagonalInverse{2,1}=diagonalInverse{1,2};%invertBlockMatrix(  subtractBlockDiagMat(B11, multiplyBlockDiagMat( B12, multiplyBlockDiagMat(A12inv, B21,decFac-1 ), decFac-1 ))  ,decFac-1)  ;
    diagonalInverse{2,2}=minusBlockDiagMat(temp1,decFac/2);




end;