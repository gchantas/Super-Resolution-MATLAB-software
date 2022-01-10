function I1=imdownsampleblur(I,d,offset,H)

I1 = real(ifft2(fft2(I).*H));

I1=downsample(I1,d,offset(1));
I1=downsample(I1',d,offset(2));
I1=I1';