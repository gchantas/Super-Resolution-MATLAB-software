function diagonalInverse=invertPseudoBlockMatrix(A,decFac,tol)

%This to work should all elementary blocks matrices be -----> diagonal

if decFac==1
    diagonalInverse=gpuArray2(zeros(size(A)));
    diagonalInverse(A<tol)=0;
    diagonalInverse(A>=tol)=gpuArray2(1./(A(A>=tol)));
   % diagonalInverse=gpuArray(conj(A)./(abs(A).^2));
else

    A11=A{1,1};
    A12=A{1,2};
    A21=A{2,1};
    A22=A{2,2};

    decFac=decFac/2;
    A11inv=invertPseudoBlockMatrix(A11,decFac,tol);
    %A12inv=invertBlockMatrix(B12,decFac-1);
    %A21inv=A12inv;%invertBlockMatrix(B21,decFac-1);
    A22inv=invertPseudoBlockMatrix(A22,decFac,tol);

    temp2=invertPseudoBlockMatrix(  subtractBlockDiagMat( multiplyBlockDiagMat(multiplyBlockDiagMat( A12, A22inv, decFac ), A21, decFac ), A11,decFac ) ,    decFac,tol);
    
    diagonalInverse{1,1}= minusBlockDiagMat(temp2,decFac);%invertBlockMatrix(  subtractBlockDiagMat(A11, multiplyBlockDiagMat(multiplyBlockDiagMat( A12, A22inv, decFac ), A21, decFac ) ,decFac )  ,decFac);%
    temp1= invertPseudoBlockMatrix( subtractBlockDiagMat( multiplyBlockDiagMat(   A21,  multiplyBlockDiagMat( A11inv , A12 , decFac), decFac) , A22, decFac), decFac,tol);
    diagonalInverse{1,2}=  multiplyBlockDiagMat(  A11inv,   multiplyBlockDiagMat(  A12,  temp1,   decFac),   decFac);
     
 %  diagonalInverse{2,1}=  multiplyBlockDiagMat(multiplyBlockDiagMat( temp1, A21, decFac )  ,  A11inv,  decFac   );
   diagonalInverse{2,1}=  multiplyBlockDiagMat( A22inv, multiplyBlockDiagMat(A21,   temp2,    decFac )  ,  decFac   );
    diagonalInverse{2,2}=minusBlockDiagMat(temp1,decFac);


end;