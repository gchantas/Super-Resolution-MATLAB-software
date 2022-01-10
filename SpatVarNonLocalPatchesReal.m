clear all
randn('seed',0);


experimentN=3;

BSNR=0;
switch experimentN
    case 1 
        blurtype=1;
        sigma=sqrt(2)/255; 
    case 2 
        blurtype=1;
        sigma=sqrt(8)/255;
    case 3 
        blurtype=4;
        sigma=sqrt(0.148)/255;
    case 4         
        blurtype=3;
        sigma=sqrt(49)/255;
    case 5         
        blurtype=4;
        siggauss=1.6;
        sigma=sqrt(4)/255;
    case 6
        blurtype=4;
        siggauss=0.4;
        sigma=sqrt(64)/255;
end


ISNRList=[];

%name ='barbara.png';
%name ='Lena512.png';
%name ='boat.png';
name = '0001.png';
decFactor=2;

f_= im2double(rgb2gray(imread(name))) ;

[Nx, Ny]=size(f_);
Nx=1024;
Ny=1024;
f=f_(100:Nx+99,100:Nx+99);


%g(1:min(Nx,Nx1),1:min(Ny,Ny1))=f1(500:499+min(Nx,Nx1),1:min(Ny,Ny1),1);

 N=Nx;
 NN = Nx*Ny;

 hh=zeros(Nx,Ny);


%     sigall=1;
%     
%     %Gaussian Blurring Function
%     sigblurx=sigall;
%     sigblury=sigall;
%     for i=1:N
%         for j=1:N
%             hh(i,j)=exp(-(i-N/2-1).^2/sigblurx).*exp(-(j-N/2-1).^2/sigblury);
%             
%         end
%     end
% %  for i=-7:7
% %         for j=-7:7
% %             hh(i+N/2+1,j+N/2+1)=1/(1+i^2+j^2);
% %         end
% %     end
%     hh=fftshift(hh)/sum(sum(hh));
%     V=fft2(hh);
        
%%%Generate BLURING to produce observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hh=zeros(Nx,Ny);

if blurtype==1
        disp('Radial blur used');
    %Generate radial blur
    for i=-7:7
        for j=-7:7
            if (i^2+j^2<30)
                hh(i+Nx/2+1,j+Ny/2+1)=1;%s/(1+i^2+j^2);
            end
        end
    end

elseif blurtype==2
    disp('Uniform blur used');
    % %Generate Uniform blur
    R=1;
    hh(Nx/2-R+1:Nx/2+R+1,Ny/2-R+1:Ny/2+R+1)=1;
    %%%%%%%%%%%%%%%%%%%%%

elseif blurtype==3
    
    disp('Pyramidal blur used');
    %Generate pyramidal blur
    hf1=zeros(Nx,Ny);
    hf2=zeros(Nx,Ny);
    
    hf1(1,1)=6/16;
    hf1(1,2)=4/16;
    hf1(1,3)=1/16;
    hf1(1,Ny)=4/16;
    hf1(1,Ny-1)=1/16;
    
    hf2(1,1)=6/16;
    hf2(2,1)=4/16;
    hf2(3,1)=1/16;
    hf2(Nx,1)=4/16;
    hf2(Nx-1,1)=1/16;
    
    hf=real(ifft2(fft2(hf1).*fft2(hf2)));
    mask=zeros(Nx,Ny);
    mask(1:3,1:3)=1;
    mask(Nx-1:Nx,Ny-1:Ny)=1;
    mask(1:3,Ny-1:Ny)=1;
    mask(Nx-1:Nx,1:3)=1;
    hh=hf.*mask; 
    
elseif blurtype==4
      sigall=1.1;
    disp('Gaussian blur used');
    %Gaussian Blurring Function
    sigblurx=sigall;
    sigblury=sigall;
    for i=1:Nx
        for j=1:Nx
            hh(i,j)=exp(-(i-Nx/2-1).^2/sigblurx).*exp(-(j-Ny/2-1).^2/sigblury);
            
        end
    end
%     
%  hh=zeros(N);
%  v=fspecial('gaussian', 35,siggauss);
%  hh(N/2-12+1:N/2+12+1,N/2-12+1:N/2+12+1)=v;

  
else disp('Wrong blurtype argument');return;
end

hh=hh/sum(sum(hh));
if blurtype~=3
    hh=fftshift(hh);
end

H=fft2(hh);

if(BSNR~=0)
    sigma=sqrt( norm( real(ifft2(H.*fft2(f))-sum(sum(f))/NN),'fro')^2/(NN*10^(BSNR/10)))
end

%PSNR=22;
%sigma=sqrt(1/10^(PSNR/10))



%%%%%%%%%%%%%%%%%%%%
%Generate Q-regularization operator

q=zeros(Nx,Ny);
q(1,1)=-4;
q(1,2)=1;
q(2,1)=1;
q(N,1)=1;
q(1,N)=1;


Q=fft2(q);


%%

% circmask=zeros(N,N);
% for k=1:N
%     for l=1:N
%         if  (k-N/2-1)^2+(l-N/2-1)^2+(k-N/2-1)*(l-N/2-1)<4
%             circmask(k,l)=1;
%         end
%     end
% end
%     circmask=fftshift(circmask);
%     


%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generate blurred data
ssigma=sigma^2;
%n=(sigma)*randn(N/decFactor);
%g=imresize(real(ifft2(H.*fft2(f))),1.0/decFactor)+n;
name = '0001x2.png';

g_ = im2double(rgb2gray(imread(name)));

g=g_(100/2:min([Nx Ny])/2+100/2-1,100/2:min([Nx Ny])/2+100/2-1);
%%%%%%%%%%%%%%%%%%%%%%
%PSNRest=10*log10( 1/(ssigma*255) )


expNum=100;
coord=zeros(expNum,2);


optcoord=zeros(Nx,Ny);
indx=1;
%  qa=zeros(N);
%  qa(1,1)=1;
%round(2*log(a+1))
for j1=1:3:15
    for k1=-j1:round(sqrt(j1)):j1
        for l1=-j1:round(sqrt(j1)):j1
 %  for k1=-j1:j1
  %    for l1=-j1:j1
            %(k1>0 | (k1==0 & l1>0))&
          if (k1>0 | (k1==0 & l1<0))&(abs(k1)==j1 | abs(l1)==j1)%(k1~=0 | l1~=0)%sqrt(k1^2+l1^2)<6&

               % l=l1
               k1+N/2+1
               l1+N/2+1
                coord(indx,1)=mod(l1,Nx);
                coord(indx,2)=mod(k1,Ny);
                optcoord(l1+Nx/2+1,k1+Ny/2+1)=1;
                
                indx=indx+1;
         end
           
         if k1==0 && l1==0
                   optcoord(l1+Nx/2+1,k1+Ny/2+1)=0.5;
         end

        end
    end
end

  expNum=indx-1
  hh=zeros(Nx,Ny);

  sigall=0.15;
  %Gaussian Blurring Function
  sigblurx=sigall;
  sigblury=sigall;

  for i=1:N
      for j=1:N
          hh(i,j)=exp(-(i-Nx/2-1).^2*sigblurx).*exp(-(j-Ny/2-1).^2*sigblury);
      end
  end

  S{1}=17;
  hh=S{1}*(hh)/sum(sum(hh));
  HQ{1}=fft2(fftshift(hh));

  g1=hh(Nx/2:Nx/2+2,Ny/2:Ny/2+2);
  sumg=sum(sum(g1));
  g1=g1/sumg;
  C=zeros(9);
  for k=1:9, C(k,k)=g1(k); end
  C1=C-g1(:)*g1(:)';
  v=eig(C1);
  v(1)=1;
  logdetC{1}=prod(sumg*v)/sumg;
  
  hh=zeros(Nx,Ny);
  sigall=0.4;
  %Gaussian Blurring Function
  sigblurx=sigall;
  sigblury=sigall;

  for i=1:N
      for j=1:N
          hh(i,j)=exp(-(i-Nx/2-1).^2*sigblurx).*exp(-(j-Ny/2-1).^2*sigblury);
      end
  end

  S{2}=23;
  hh=S{2}*(hh)/sum(sum(hh));
  HQ{2}=fft2(fftshift(hh));
  g1=hh(Nx/2:Nx/2+2,Ny/2:Ny/2+2);
  sumg=sum(sum(g1));
  g1=g1/sumg;
  C=zeros(9);
  for k=1:9, C(k,k)=g1(k); end
  C1=C-g1(:)*g1(:)';
  v=eig(C1);
  v(1)=1;
  logdetC{2}=prod(sumg*v)/sumg;
 
  
  sigall=2;
  %Gaussian Blurring Function
  sigblurx=sigall;
  sigblury=sigall;
  hh=zeros(N);
  for i=1:N
      for j=1:N
          hh(i,j)=exp(-(i-Nx/2-1).^2*sigblurx).*exp(-(j-Ny/2-1).^2*sigblury);
      end
  end

  S{3}=13;
  hh=S{3}*(hh)/sum(sum(hh));
  HQ{3}=fft2(fftshift(hh));

  g1=hh(Nx/2:Nx/2+2,Ny/2:Ny/2+2);
  sumg=sum(sum(g1));
   g1=g1/sumg;
  C=zeros(9);
  for k=1:9, C(k,k)=g1(k); end
  C1=C-g1(:)*g1(:)';
  v=eig(C1);
  v(1)=1;
  logdetC{3}=prod(sumg*v)/sumg;

  
  
  
  sigall=0.2;
  %Gaussian Blurring Function
  sigblurx=sigall;
  sigblury=sigall;
  hh=zeros(N);
  for i=1:N
      for j=1:N
          hh(i,j)=exp(-(i-Nx/2-1).^2*sigblurx).*exp(-(j-Ny/2-1).^2*sigblury);
      end
  end

  S{4}=30;
  hh=S{4}*(hh)/sum(sum(hh));
  HQ{4}=fft2(hh);
  
  g1=hh(Nx/2:Nx/2+2,Ny/2:Ny/2+2);
  sumg=sum(sum(g1));
    g1=g1/sumg;
  C=zeros(9);
  for k=1:9, C(k,k)=g1(k); end
  C1=C-g1(:)*g1(:)';
  v=eig(C1);
  v(1)=1;
  logdetC{4}=prod(sumg*v)/sumg;
  
  kerNum=1;
  
    
%[fhat,a,ssigma]=EM(single(g),single(real(ifft2(H))),single(real(ifft2(Q))),N,N,100);
%[a,ssigma,fhat]=stat_rest(1,1,g,Q,H,N);
% a=50000;
%ssigma

 Hf = abs(H).^2;

 Qf =ones(Nx,Ny);%abs(Q).^2;

%Mfg=conj(H).*fft2(g)./(Hf+a*ssigma*Qf/10);
%fhat=real(ifft2(Mfg));

%isnrstat=10*log10(norm(f-g,'fro')^2/norm(f-fhat,'fro')^2)

%PSNRtrue=10*log10( NN/norm(g-f,'fro')^2 )



%Qp=zeros(N);
   
   for k=1:kerNum
       K{k}=ones(Nx,Ny)/kerNum;
   end
   
   
%==================Start Iterations==========================

%%
x1=imresize( g,  decFactor  );


%  hhz=zeros(N);
%  v=fspecial('gaussian', 25,0.5);
%  hhz(N/2-12+1:N/2+12+1,N/2-12+1:N/2+12+1)=v;
%  HHZ=fft2(fftshift(hhz));

DFTShift=ones(N,N);

sigma_method=1;

P=1;
N_low=Nx/decFactor;
Lall=zeros(N_low,N_low,P);
Lall(:,:,1)=g;
dec = zeros(P,2); %This should be zero

psigma = sigma*ones(P,1);
dall   = zeros(P-1,2);
%EM using stationary prior
x=x1;

a=1;
dfac=decFactor;

for iter=1:50

    DTemp(:,:,1)=DFTShift(:,:,1)/(psigma(1));
    Y1=fft2( imupsampleblur(Lall(:,:,1), decFactor, dec(1,1:2), H(:,:,1))).*DTemp(:,:,1)/psigma(1);
    for k=2:P
      [dummy, DFTShift(:,:,k)]=dft_shift(zeros(N,N),dall(k-1,:),N);
      DTemp(:,:,k)=DFTShift(:,:,k)/(psigma(k));
      Y1=Y1+fft2(imupsampleblur(Lall(:,:,k), dfac, dec(k,1:2),H(:,:,k))).*conj(DTemp(:,:,k)/psigma(k));
    end
    
    xprev=x1;
    
    [sum1,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I] = traceCom(ones(N,N),DTemp.*H,abs(Q).^2,a,1,N,P,abs(Q).^2);
    %Analytical inversion
    X = invertY(Y1,N,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I);
    x1 = real(ifft2( X ));
    a = (NN-1)/(sum(sum(real(ifft2( Q.*X )).^2))+sum1)

    %One sigma for each image
    if sigma_method==1
     
        s1 = real(HCfgCom(ones(N,N),DFTShift.*H,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I,N,P));
        s2=0;
 
        for k=1:P
            s2 = s2 + norm(Lall(:,:,k)-imdownsampleblur(x1,dfac,dec(k,1:2),H(:,:,k).*DFTShift(:,:,k)),'fro')^2;
        end
        
        ssigma=(s1+s2)/(P*NN/4)
        psigma(:)=sqrt(ssigma);
    end
    
    %Multiple sigma
    if sigma_method==2


        for k=1:P

            TEMP=ones(N,N,P);
            for m=1:P
                TEMP(:,:,m)=H(:,:,k).*DTemp(:,:,k)/sqrt(P);
            end
            
            s1 = real(HCfgCom(ones(N,N),TEMP,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I,N,P))
            s2 = norm(  Lall(:,:,k)-imdownsampleblur(x1,dfac,dec(k,1:2),H(:,:,k).*DFTShift(:,:,k)),'fro'   )^2
            psigma(k)=sqrt((s1+s2)/(NN/4));

        end
        psigma.^2
    end


    for k=2:P
        [dall(k-1,:) DFTShift(:,:,k)] = shift_Newton(fft2(Lall(:,:,k)),fft2(x1),dall(k-1,:),H(:,:,k),N);
    end

end


g2=imresize(g,decFactor);
ISNR_Bayes=10*log10(norm(f(20:Nx-20,20:Ny-20)-g2(20:Nx-20,20:Ny-20),'fro')^2/norm(f(20:Nx-20,20:Ny-20)-x1(20:Nx-20,20:Ny-20),'fro')^2)


for exper=1:expNum
    nu=1;
    c2{exper}=100;
    rw{2*exper-1}(1)=coord(exper,2);
    rw{2*exper-1}(2)=coord(exper,1);
    Z{exper}=ones(Nx,Ny)/expNum;
    pof{exper}=1/expNum;
   Qp=zeros(Nx,Ny);
end;

for iter=1:2

    ZALL=zeros(Nx,Ny);
    
%     SA2=zeros(N);
% 
%     for exper=1:expNum
% 
%         SA=zeros(N);
%         qa=zeros(N);
%         qa(1,1)=1;
%         qa(coord(exper,2)+1,coord(exper,1)+1)=-1;
%         q1=real(ifft2(fft2(qa))).^2;
% 
% %         if iter>1
% %             KALL=zeros(N);
% %             for k=1:kerNum
% %                 E =  real(ifft2( fft2(  (x-circshift(x,rw{2*exper-1}) ) .^2 + real(ifft2(fft2(q1).*(fft2(Qp)))) ) .* HQ{k}  ));
% %                % E =  real(ifft2( fft2(  (x-circshift(x,rw{2*exper-1}) ) .^2 + Qp ) .* HQ{k}  ));
% %                 K{k}=exp(Z{exper}.*(  -0.5*S{k}*log(2*pi*exp(1)) - S{k}*log(2*pi/c2{exper}) + 0.5*S{k}*log(A{exper}) - 0.5*c2{exper}*A{exper}.*E ));
% %                % K{k}=exp(Z{exper}.*( -0.5*S{k}*log(2*pi/c2{exper}) + 0.5*S{k}*log(A{exper}) - 0.5*c2{exper}*A{exper}.*E ));
% %                 KALL=KALL+K{k};
% %             end
% %             for k=1:kerNum
% %                 K{k}=K{k}./KALL;
% %             end
% % %             TEMP=K{1};
% % %             K{1}=K{2};
% % %             K{2}=TEMP;
% %         end
% 
%         E=zeros(N);
%         
%         for k=1:kerNum
%             E =  real(ifft2( fft2(  (x-circshift(x,rw{2*exper-1}) ) .^2 + real(ifft2(fft2(q1).*(fft2(Qp)))) ) .* HQ{k}  ));
%             SA2=SA2+Z{exper}.*K{k}*S{k};
%             SA=SA+K{k}*S{k};
%         end
% 
%         J=(SA/2+nu/2);
%         A{exper} = J./((c2{exper}/2)*E+nu/2);
%         Z{exper}=(nu+c2{exper}*E).^(-J);
%         %Z{exper} = pof{exper}*exp( - 0.5*SA * log(2*pi/c2{exper})+J.*log(A{exper}) - A{exper} .* (c2{exper}*E/2+nu/2) );
%         ZALL=ZALL+Z{exper};
%     end
% 
%    % SA2=zeros(N);
% 
%     for exper=1:expNum
%         Z{exper}=Z{exper}./ZALL;
%        % pof{exper}=sum(sum(Z{exper}))/NN;
%  %
% %         if iter>1
% %             KALL=zeros(N);
% %             for k=1:kerNum
% %                 E =  real(ifft2( fft2(  (x-circshift(x,rw{2*exper-1}) ) .^2 + real(ifft2(fft2(q1).*(fft2(Qp)))) ) .* HQ{k}  ));
% %                 %K{k}=exp(Z{exper}.*(  -0.5*S{k}*log(2*pi*exp(1)) - S{k}*log(2*pi/c2{exper}) + 0.5*S{k}*log(A{exper}) - 0.5*c2{exper}*A{exper}.*E ));
% %                 K{k}=exp(Z{exper}.*( -0.5*S{k}*log(2*pi/c2{exper}) + 0.5*S{k}*log(A{exper}) - 0.5*c2{exper}*A{exper}.*E ));
% %                 KALL=KALL+K{k};
% %             end
% %             for k=1:kerNum
% %                 K{k}=K{k}./KALL;
% %             end
% %         end
% % 
% %         for k=1:kerNum
% %             SA2=SA2+Z{exper}.*K{k}*S{k};
% %         end
% 
%        % Z{exper}(find(Z{exper}<0.05))=0.0;
%         %B{exper}=B{exper}./ZALL;
% 
%         B{exper}=zeros(N);
%         for k=1:kerNum
%             B{exper} = B{exper} + c2{exper}*real(ifft2( fft2( K{k} .* Z{exper} .* A{exper} ) .* HQ{k} ))./SA;
%         end
%         B{exper} = single( B{exper} );
%     end

     for exper=1:expNum

        J=(S{kerNum}/2+nu/2);

        qa=zeros(Nx,Ny);
        qa(1,1)=1;
        qa(coord(exper,2)+1,coord(exper,1)+1)=-1;


        q1=real(ifft2(fft2(qa))).^2;
        E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
        E{exper} = real(ifft2( fft2( ( x-circshift(x, rw{2*exper-1})  ).^2  + real(ifft2(fft2(q1).*(fft2(Qp))))) .* HQ{kerNum}  ));
        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
        A{exper} = (nu+S{kerNum})./(c2{exper}*E{exper}+nu);
        Z{exper} = exp(  J*log(A{exper})-A{exper}.*(c2{exper}*E{exper}/2+nu/2)   );
        %Z{exper}=pof{exper}.*(nu+c2{exper}*E).^(-J);
        %B{exper}=A{exper};
        ZALL=ZALL+Z{exper};

    end

    for exper=1:expNum
        Z{exper}=Z{exper}./ZALL;
       % Z{exper}(find(Z{exper}<0.05))=0.0;
        B{exper}=A{exper}.*Z{exper};
        B{exper}= ssigma*c2{exper}*real(ifft2( fft2( B{exper}).* HQ{kerNum} ))/S{kerNum};
        rw{2*exper}=B{exper};
    end
    
   % A0=zeros(N);
%     for exper=1:expNum
%         
%         qa=zeros(N);
%         qa(1,1)=1;
%                 
%         qa(coord(exper,2)+1,coord(exper,1)+1)=-1;
%         q1=real(ifft2(fft2(qa))).^2;
%         %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + real(ifft2(fft2(q1).*(fft2(Qp))))) .*( HQ)  )));
%        
%         A = (nu+S)./(c2{exper}*E+nu);
%         
%          %A0 = A0 + Z{exper}+10000;
%         %Z{exper}=ones(N);
%          c2{exper}=S*sum(sum(Z{exper}));
%          c2{exper}=c2{exper}/sum(sum((Z{exper}.*A).*E));
%         %c1{exper}=bisection(0.00000001,100,B{exper}./ZALL,Z{exper},c1{exper},S,N);
%     end

    Ainv=single( real(ifft2(zeros(Nx,Ny))) );
    xprev=x;
    
    Qp=sum(sum( real(ifft2( H )).*real(ifft2( conj(H) )) ))/ssigma;

    for exper=1:expNum

        %pof{exper}=exp( psi(Z{exper}+10000) -psi(A0) );%(Z{exper}+0.0001)./ A0;
        qa=zeros(Nx,Ny);
        qa(1,1)=1;

        qa(coord(exper,2)+1,coord(exper,1)+1)=-1;
        q1=real(ifft2(fft2(qa))).^2;
        Qp=Qp+real(ifft2(fft2(q1).*(fft2(B{exper}))));
    end

    Qp=0./Qp;
    
    %  Qp=sum(sum( abs(H).^2 ))/(ssigma);
    %
    %     for exper=1:expNum
    %
    %         qa=zeros(N);
    %         qa(1,1)=1;
    %
    %         qa(coord(exper,2)+1,coord(exper,1)+1)=-1;
    %
    %         Qp=Qp+sum(sum(abs(fft2(qa)).^2))*sum(sum(B{exper}));
    %     end
    %
    %     Qp=Qp/NN;
    %%%%%%Image Estimation

   % tic
%    x  = cudaSymLanczosBlockMatch(single(zeros(N)),single(real(ifft2(H))/(sqrt(ssigma))),expNum,coord',B,single(ones(N)),single( real(ifft2(conj(H).*fft2(g)))/(ssigma)) ,1,N,N,1,0.01,150);
    % b=real(ifft2(conj(H).*fft2(g)))/(ssigma);
     %iw=D.^2/ssigma;
    
    % x=reshape(cgLanczos(@Amat,zeros(NN,1), b(:), iw, rw, 1, 0, 550, 0.0001),N1,N2);

    %toc
%     iter
%     if iter == 3 || iter==6
%         x=x+0.01*randn(N);
%     end
    itnum=0;
%         for exper=1:expNum
%             Qp{exper}=zeros(N);
%         end

    if itnum>0

        for iter1=1:itnum
            randn('state',iter1);
            d=randn(Nx,Ny);
            
            %e  = cudaSymLanczosBlockMatch(single(zeros(N)),single(real(ifft2(H))/sqrt(ssigma)),expNum,coord',B,single(ones(N)),single(d),ssigma,N,N,R,0.0001,150) ;
            iw=abs(H).^2/ssigma;
            b=d;
            tic;e=reshape(cgLanczos(@Amat,zeros(NN,1), b(:),iw,rw, 1, 0, 1050, 0.000001),N,N);toc

            for exper=1:expNum
                Qp{exper} = Qp{exper} +  (e-circshift(e,rw{2*exper-1})) .*  (d-circshift(d,rw{2*exper-1}))/itnum;
            end
        end

    end

    iw{1}=H;
    iw{2}=decFactor;
    g1=zeros(Nx,Ny);%imresize(  g,   decFactor);
    %g1(2:decFactor:Nx, 2:decFactor:Ny)=g;
    
    
    g1 = upsample(g,2,0);
    g1 = upsample(g1',2,0);
    g1=g1';

    b = real(ifft2( conj(H) .* fft2( g1 ) ));

    tic;x=reshape(cgLanczos(@Amat,zeros(NN,1), b(:), iw, rw, 1, 0,1000, 0.000000001), Nx, Ny);toc
    norm(x-xprev,'fro')/NN

    %ISNR_Bayes=10*log10(norm(f-imresize(g,decFactor),'fro')^2/norm(f-x,'fro')^2)

    ISNR_Bayes=10*log10(norm(f(20:Nx-20,20:Ny-20)-g2(20:Nx-20,20:Ny-20),'fro')^2/norm(f(20:Nx-20,20:Ny-20)-x(20:Nx-20,20:Ny-20),'fro')^2)

  % options = optimoptions(@fminunc,'Algorithm','quasi-newton');
  % sigall = bisectionSigma(sigall,Z{1},A{1},E1{1},N);
  % Gaussian Blurring Function
  %   sigblurx=sigall;
  %   sigblury=sigall;
  %   hh=zeros(N);
  %   for i=1:N
  %       for j=1:N
%           hh(i,j)=exp(-(i-N/2-1).^2*sigblurx).*exp(-(j-N/2-1).^2*sigblury);
%       end
%   end
%   S{1}=7;
%   hh=S{1}*(hh)/sum(sum(hh));
%   HQ{1}=fft2(fftshift(hh));
  
    if iter==100
        U=zeros(N);
        clear Z;
        clear A;
        
        ZALL=zeros(N);
        for exper=1:expNum
            E = (real(ifft2( fft2( ( f-circshift(f,rw{2*exper-1})  ).^2  ) .* (HQ{exper} ) )));
            A = (nu+S{exper})./(nu*E+nu);
            Z{exper}=exp(J*(psi(J)-log(J)+log(A))-A.*(E/2+nu/2));
           % Z{exper}=(gamma(J)/gamma(c1{exper}/2))*(c2{exper}/c1{exper})^(S/2)*(1+c2{exper}*E/c1{exper}).^(-J);
            ZALL=ZALL+Z{exper};
        end
%         for exper=1:expNum
%             
%             E =  ( f-circshift(f,rw{2*exper-1}) ).^2   ;
%             A = (nu+S{exper})./(real(ifft2( fft2( E ).*(HQ{exper}))) +nu);
%             Z{exper}=Z{exper}./ZALL;
%             %
%             %              %Z{exper}=(gamma(J)/gamma(c1{exper}/2))*(c2{exper}/c1{exper})^(S/2)*(1+c2{exper}*E/c1{exper}).^(-J);
%             %
%             %              % Z(find(real(ifft2(fft2(Z.*HQ{exper})))<0.001))=0.000000001;
%             %
%             %              B{exper}=Z{exper}.*A;
%         
%             %U = real(ifft2( conj(fft2( (fftshift( (A.*Z{exper}./ZALL))))) .* fft2(E) )) ;
% 
%             R=2;
% 
%            
%             m=zeros(2*R+1,2*R+1);
%             U=zeros(2*R+1);
%             h=zeros(2*R+1);
%             for i=R+1:N-R-1
%                 for j=R+1:N-R-1
%                     if Z{exper}(i,j)>0.1
%                     U=0.5*E(i-R:i+R,j-R:j+R)*A(i,j);
%                     U=U-min(min(U));
%                     h=(exp(-U))/(sum(sum( exp(-U) )));
%                     
%                     m=m+Z{exper}(i,j)*h/sum(sum((h)));
%                    
%                     
% 
%                   %  figure(1),imagesc(h),colormap('gray');
%                     end
%                 end
%             end
%             
%             
%             
% 
%             h2=zeros(N);
%             h2(N/2-R+1:N/2+R+1,N/2-R+1:N/2+R+1) = m ;
%             h2=fftshift(h2);
%             h2=S{exper}*h2/sum(sum(h2));
%             
%             HQ{exper}=(fft2(h2));
%             
%         end
        
    end

 end




