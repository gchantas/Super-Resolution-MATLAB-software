function [sum1,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I]=traceCom(H,DFTShift,Q,a,ssigma,Nx,Ny,P,M)
NN=Nx*Ny;
k=1:Ny/2;
k1=1:Ny/2;
sum1=0;
A11I = []; A12I = []; B11I = []; B12I = [];
A21I = []; A22I = []; B21I = []; B22I = [];
A1I = []; A2I = []; B1I = []; B2I = [];
C1I = []; C2I = []; D1I = []; D2I = [];
dxy=2;

S = zeros(P,NN);

for m=1:P
    Sm=DFTShift(:,:,m)/dxy;
    S(m,:) = H(:).*Sm(:);
end

for n=1:Nx/2

    HS(k1)=a*ssigma.*Q(k);
    for m=1:P
        HS(k1) = HS(k1)+abs(S(m,k)).^2;
    end

    HS(k1+Ny/2)=a*ssigma.*Q(k+Ny/2);
    for m=1:P
        HS(k1+Ny/2) = HS(k1+Ny/2)+abs(S(m,k+Ny/2)).^2;
    end


    H1=0;
    for m=1:P
        H1 = H1 + conj(S(m,k)).*S(m,Ny/2+k);
    end
%H1 = conj(S1(k)).*S1(N/2+k) + conj(S2(k)).*S2(N/2+k) + conj(S3(k)).*S3(N/2+k) + conj(S4(k)).*S4(N/2+k);

TEMP = (HS(Ny/2+k1).*HS(k1)-abs(H1).^2);
Ai1 = HS(Ny/2+k1)./TEMP;
Ai2 = HS(k1)./TEMP;

Bi1 = -H1./TEMP;
Bi2 = conj(Bi1);
clear TEMP;

   HS2(k1)=a*ssigma.*Q(k+NN/2);
    for m=1:P
        HS2(k1) = HS2(k1) + abs(  S(m,  NN/2+k)   ).^2;
    end
    
    HS2(k1+Ny/2) = a*ssigma.*Q(NN/2+Ny/2+k);
    for m=1:P
        HS2(  k1  +   Ny/2   ) = HS2(k1+Ny/2)  +  abs(  S(m,NN/2+Ny/2+k)   ).^2;
    end


    H12=0;
    for m=1:P
        H12 = H12 + conj(S(  m,   NN/2+k  )).*S(  m,  NN/2+Ny/2+k   );
    end
    
    
%HS2(k1) = abs(S1(NN/2+k)).^2 + abs(S2(NN/2+k)).^2 + abs(S3(NN/2+k)).^2 + abs(S4(NN/2+k)).^2+a*ssigma.*Q(NN/2+k);
%HS2(N/2+k1) = abs(S1(NN/2+N/2+k)).^2 + abs(S2(NN/2+N/2+k)).^2 + abs(S3(NN/2+N/2+k)).^2 + abs(S4(NN/2+N/2+k)).^2+a*ssigma.*Q(NN/2+N/2+k);
%H12 = conj(S1(NN/2+k)).*S1(NN/2+N/2+k) + conj(S2(NN/2+k)).*S2(NN/2+N/2+k) + conj(S3(NN/2+k)).*S3(NN/2+N/2+k) + conj(S4(NN/2+k)).*S4(NN/2+N/2+k);

TEMP = (HS2(Ny/2+k1).*HS2(k1)-abs(H12).^2);
Ai12 = HS2(Ny/2+k1)./TEMP;
Ai22 = HS2(k1)./TEMP;

Bi12 = -H12./TEMP;
Bi22 = conj(Bi12);

C1=0;
for m=1:P
    C1 = C1 + conj(S(m,k)).*S(  m, NN/2+k   );
end
C2=0;
for m=1:P
    C2 = C2 + conj(S(m,k+Ny/2)).*S(   m,  NN/2+k+Ny/2   );
end

D1=0;
for m=1:P
    D1 = D1 + conj(S(m,k)).*S(m,  NN/2+Ny/2+k   );
end

D2=0;
for m=1:P
    D2 = D2 + conj(S(m,k+Ny/2)).*S(m, NN/2+k);
end

%
C211 = Ai1.*abs(C1).^2+2*real(C1.*conj(D2).*Bi2)+Ai2.*abs(D2).^2;
C222 = Ai1.*abs(D1).^2+2*real(C2.*conj(D1).*Bi1)+Ai2.*abs(C2).^2;
C212 = Ai1.*conj(C1).*D1+Bi2.*conj(D2).*D1+Bi1.*C2.*conj(C1)+conj(D2).*Ai2.*C2;
C221 = conj(C212);
%

TEMP = (HS2(Ny/2+k1)-C222).*(HS2(k1)-C211) - (H12-C212).*(conj(H12)-C221);
AINV21 = (HS2(Ny/2+k1)-C222)./ TEMP;
AINV22 = (HS2(k1)-C211)./ TEMP;
BINV21 = -(H12-C212)./ TEMP;
BINV22 = conj(BINV21);
%

C1=conj(C1);
C2=conj(C2);
DT=D2;
D2=conj(D1);
D1=conj(DT);

C211 = Ai12.*abs(C1).^2+2*real(C1.*conj(D2).*Bi22)+Ai22.*abs(D2).^2;
C222 = Ai12.*abs(D1).^2+2*real(C2.*conj(D1).*Bi12)+Ai22.*abs(C2).^2;
C212 = Ai12.*conj(C1).*D1+Bi22.*conj(D2).*D1+Bi12.*C2.*conj(C1)+conj(D2).*Ai22.*C2;
C221 = conj(C212);

TEMP = (HS(Ny/2+k1)-C222).*(HS(k1)-C211) - (H1-C212).*(conj(H1)-C221);
AINV11 = (HS(Ny/2+k1)-C222)./ TEMP;
AINV12 = (HS(k1)-C211)./ TEMP;
BINV11 = -(H1-C212)./ TEMP;
BINV12 = conj(BINV11);

C111 = HS2(k1);
C122 = HS2(k1+Ny/2);
C112 = H12;
C121 = conj(C112);

C1=conj(C1);
C2=conj(C2);
DT=D2;
D2=conj(D1);
D1=conj(DT);

% TEMP = (HS2(N/2+k1).*HS2(k1)-abs(H12).^2);
% Ai12 = HS2(N/2+k1)./TEMP;
% Ai22 = HS2(k1)./TEMP;
% 
% Bi12 = -H12./TEMP;
% Bi22 = conj(Bi12);

K1 = AINV11.*C1+BINV11.*D2;
K2 = AINV11.*D1+BINV11.*C2;

%C1K = (K1.*C122./C121-K2)./(C112-C111.*C122./C121);
C1K = -K1.*Ai12-K2.*Bi22;
D1K = (-K1-C1K.*C111)./C121;

K1 = BINV12.*C1+AINV12.*D2;
K2 = BINV12.*D1+AINV12.*C2;

%D2K = (K1.*C122./C121-K2)./(C112-C111.*C122./C121);
D2K = -K1.*Ai12-K2.*Bi22;
C2K = (-K1-D2K.*C111)./C121;

sum1=sum1+sum(  AINV11.*M(k)   );
sum1=sum1+sum(  AINV12.*M(k+Ny/2)   );

sum1=sum1+sum( AINV21.*M(k+NN/2) );
sum1=sum1+sum( AINV22.*M(k+NN/2+Ny/2) );


A11I=[A11I AINV11];
A12I=[A12I AINV12];
B11I=[B11I BINV11];
B12I=[B12I BINV12];

A21I=[A21I AINV21];
A22I=[A22I AINV22];
B21I=[B21I BINV21];
B22I=[B22I BINV22];


C1I=[C1I C1K];C2I=[C2I C2K];
D1I=[D1I D1K];D2I=[D2I D2K];

k=k+Ny;
end;

A1I = [A11I A21I];
A2I = [A12I A22I];

B1I = [B11I B21I];
B2I = [B12I B22I];
