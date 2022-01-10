function g = grad_gtn(n,A,nprev,c,N, M)
NN=N*M;
%mysumA=sum(sum( log(A) ))-sum(sum(A));
%lp = mysumA + NN*(psi(nprev/2+1/c)-log(nprev/2+1/c) - (psi(n/2)-1-log(n/2)));

g = (NN/2)*(psi(nprev/2+1/c)-log((nprev/2+1/c)))+sum(sum(log(A)))/2-(c*n^(c/2-1)/2^(c/2+1))*sum(sum(A))+(NN*c/4)*log(n/2)+NN*c/4-NN*psi(n/2)/2;