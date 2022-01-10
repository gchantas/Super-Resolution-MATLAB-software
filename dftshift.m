function [I1 K] = dft_shift(I,d,N)
NN=N*N;

for n=1:N, for m=1:N, 
        K(n,m) = exp(-i*2*pi*d(1)*(n-N/2-1)/N)*exp(-i*2*pi*d(2)*(m-N/2-1)/N); 
    end; end;

K=fftshift(K);

I1=real(ifft2( fft2(I).*K ));