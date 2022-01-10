function YI=invertY(Y,Nx,Ny,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I)

NN=Nx*Ny;

k=1:Ny/2;
k1=1:Ny/2;
YI=zeros(Nx,Ny);

for n=1:Nx/2
    YI(k) = Y(k).*A1I(k1) +  Y(k+Ny/2).*B1I(k1) + Y(k+NN/2).*C1I(k1) + Y(k+Ny/2+NN/2).*D1I(k1);
    YI(k+Ny/2) = Y(k).*B2I(k1) +  Y(k+Ny/2).*A2I(k1) + Y(k+NN/2).*D2I(k1) + Y(k+Ny/2+NN/2).*C2I(k1);
    YI(k+NN/2) = Y(k+NN/2).*A1I(k1+NN/4) +  Y(k+Ny/2+NN/2).*B1I(k1+NN/4) + Y(k).*conj(C1I(k1)) + Y(k+Ny/2).*conj(D2I(k1));
    YI(k+Ny/2+NN/2) = Y(k+NN/2).*B2I(k1+NN/4) +  Y(k+Ny/2+NN/2).*A2I(k1+NN/4) + Y(k).*conj(D1I(k1)) + Y(k+Ny/2).*conj(C2I(k1));
    k=k+Ny;
    k1=k1+Ny/2;
end;