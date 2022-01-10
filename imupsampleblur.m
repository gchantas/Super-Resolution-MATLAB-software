function I1=imupsample(I,d,offset,H)

I1=upsample(I,d,offset(1));
I1=upsample(I1',d,offset(2));
I1=I1';

I1 = real(ifft2(fft2(I1).*conj(H)));