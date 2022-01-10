function [alpha,ssigma,fhat,Cfg] = stat_rest(alpha,ssigma,g,Q,H1,Nx,Ny)

NN=Nx*Ny;
G=fft2(g);
Hf1 = abs(H1).^2;
Qf = abs(Q).^2;



p=100;

%==================Start Iterations==========================

for k=1:p
    Mfg=conj(H1).*(G)./(Hf1+ssigma*alpha*Qf );
    Cfg=ssigma./( Hf1+ssigma*alpha*Qf);
    s3=sum(sum( ( Cfg + (abs(Mfg).^2)/NN ).*Qf ));
    alpha=(NN-1)/s3;
    s1=sum(sum( (abs(H1).^2).*( Cfg + (abs(Mfg).^2)/NN ) ));
    s2=sum(sum( ( abs(G).^2-2*real(conj(G).*H1.*Mfg) ) ))/NN;
    ssigma=(s1+s2)/NN;
    %ISNR_Stat=10*log10(norm(f-g,'fro')^2/norm(f-real(ifft2(Mfg)),'fro')^2);
    % pause
end;

fhat = real(ifft2(Mfg));