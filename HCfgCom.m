function sum1=HCfgCom(H,DFTShift,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I,Nx,Ny,P)
NN=Nx*Ny;
dxy=2;

k=1:Ny/2;
k1=1:Ny/2;
sum1=0;

S = zeros(P,NN);

for m=1:P
    Sm=DFTShift(:,:,m)/dxy;
    size(Sm)
    size(S) 
    size(H)
    S(m,:) = H(:).*Sm(:);
end


for n=1:Nx/2

    HS(k)=0;
    for m=1:P
        HS(k) = HS(k)+abs(S(m,k)).^2;
    end

    HS(k+Ny/2)=0;
    for m=1:P
        HS(k+Ny/2) = HS(k+Ny/2)+abs(S(m,k+Ny/2)).^2;
    end

    H1=0;
    for m=1:P
        H1 = H1 + conj(S(m,k)).*S(m,Ny/2+k);
    end
% 
% HS(k) = abs(S1(k)).^2 + abs(S2(k)).^2 + abs(S3(k)).^2 + abs(S4(k)).^2;
% HS(N/2+k) = abs(S1(N/2+k)).^2 + abs(S2(N/2+k)).^2 + abs(S3(N/2+k)).^2 + abs(S4(N/2+k)).^2;
% H1 = conj(S1(k)).*S1(N/2+k) + conj(S2(k)).*S2(N/2+k) + conj(S3(k)).*S3(N/2+k) + conj(S4(k)).*S4(N/2+k);

C1=0;
for m=1:P
    C1 = C1 + conj(S(m,k)).*S(m,NN/2+k);
end
C2=0;
for m=1:P
    C2 = C2 + conj(S(m,k+Ny/2)).*S(m,NN/2+k+Ny/2);
end

D1=0;
for m=1:P
    D1 = D1 + conj(S(m,k)).*S(m,NN/2+Ny/2+k);
end

D2=0;
for m=1:P
    D2 = D2 + conj(S(m,k+Ny/2)).*S(m,NN/2+k);
end
% C1 = conj(S1(k)).*S1(NN/2+k)+conj(S2(k)).*S2(NN/2+k)+conj(S3(k)).*S3(NN/2+k)+conj(S4(k)).*S4(NN/2+k);
% C2 = conj(S1(k+N/2)).*S1(NN/2+k+N/2)+conj(S2(k+N/2)).*S2(NN/2+k+N/2)+conj(S3(k+N/2)).*S3(NN/2+k+N/2)+conj(S4(k+N/2)).*S4(NN/2+k+N/2);
% 
% D1 = conj(S1(k)).*S1(NN/2+N/2+k)+conj(S2(k)).*S2(NN/2+N/2+k)+conj(S3(k)).*S3(NN/2+N/2+k)+conj(S4(k)).*S4(NN/2+N/2+k);
% D2 = conj(S1(k+N/2)).*S1(NN/2+k)+conj(S2(k+N/2)).*S2(NN/2+k)+conj(S3(k+N/2)).*S3(NN/2+k)+conj(S4(k+N/2)).*S4(NN/2+k);

sum1 = sum1 + sum(A1I(k1).*HS(k)+B2I(k1).*H1+conj(C1I(k1)).*C1+conj(D1I(k1)).*D1) + sum(B1I(k1).*conj(H1)+A2I(k1).*HS(k+Ny/2)+conj(D2I(k1)).*D2+conj(C2I(k1)).*C2);

   HS(k)=0;
    for m=1:P
        HS(k) = HS(k)+abs(S(m,NN/2+k)).^2;
    end
    
    HS(k+Ny/2)=0;
    for m=1:P
        HS(k+Ny/2) = HS(k+Ny/2)+abs(S(m,k+NN/2+Ny/2)).^2;
    end


    H1=0;
    for m=1:P
        H1 = H1 + conj(S(m,NN/2+k)).*S(m,NN/2+Ny/2+k);
    end

% 
% HS(k) = abs(S1(NN/2+k)).^2 + abs(S2(NN/2+k)).^2 + abs(S3(NN/2+k)).^2 + abs(S4(NN/2+k)).^2;
% HS(N/2+k) = abs(S1(NN/2+N/2+k)).^2 + abs(S2(NN/2+N/2+k)).^2 + abs(S3(NN/2+N/2+k)).^2 + abs(S4(NN/2+N/2+k)).^2;
% 
% H1 = conj(S1(NN/2+k)).*S1(NN/2+N/2+k) + conj(S2(NN/2+k)).*S2(NN/2+N/2+k) + conj(S3(NN/2+k)).*S3(NN/2+N/2+k) + conj(S4(NN/2+k)).*S4(NN/2+N/2+k);
% 
sum1 = sum1 + sum( C1I(k1).*conj(C1)+D2I(k1).*conj(D2)+A1I(k1+NN/4).*HS(k)+B2I(k1+NN/4).*H1 ) + sum( D1I(k1).*conj(D1)+C2I(k1).*conj(C2)+A2I(k1+NN/4).*HS(k+Ny/2)+B1I(k1+NN/4).*conj(H1));


k=k+Ny;
k1=k1+Ny/2;

end;