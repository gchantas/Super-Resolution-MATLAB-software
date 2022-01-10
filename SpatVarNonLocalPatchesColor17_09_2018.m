%function [PSNR, ISNR, SSIMres]=SpatVarNonLocalPatchesColor17_09_2018(imNumber,myinput1)
clear all;
myinput1=1;
imNumber=801;
%randn('seed',0);

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
name = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_HR/0%d.png',imNumber);
decFactor=2;

extras=4;%expanid by $extras pixels the boundary of the image by replicating the bounds

f_= im2double(rgb2gray(imread(name))) ;
fcolor=im2double(imread(name));


[Nx, Ny,dummy]=size(f_);
  %Nx=min(1024,Nx);
  %Ny=min(1024,Ny);
%f=f_(100:Nx+99,100:Nx+99);
f=f_(1:Nx,1:Ny);

ssigma=0.00001;
%g(1:min(Nx,Nx1),1:min(Ny,Ny1))=f1(500:499+min(Nx,Nx1),1:min(Ny,Ny1),1);

 N=Nx;
 NN = Nx*Ny;


    sigall=1;
    %Gaussian Blurring Function
    sigblurx=sigall;
    sigblury=sigall;
    for i=1:N
        for j=1:N
            hh(i,j)=exp(-(i-N/2-1).^2/sigblurx).*exp(-(j-N/2-1).^2/sigblury);
            
        end
    end
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
    sigall=2.0;
    disp('Gaussian blur used');
    %Gaussian Blurring Function
    sigblurx=sigall;
    sigblury=sigall;
    for i=1:Nx
        for j=1:Nx
            hh(i,j)=exp(-(i-Nx/2-1).^2/sigblurx).*exp(-(j-Ny/2-1).^2/sigblury);
        end
    end

%  hh=zeros(N);
%  v=fspecial('gaussian', 35,siggauss);
%  hh(N/2-12+1:N/2+12+1,N/2-12+1:N/2+12+1)=v;

  
else
    disp('Wrong blurtype argument');return;
end

Nx=Nx+2*extras;
Ny=Ny+2*extras;

hh1=zeros(Nx,Ny);
sigall=0.5;
 

disp('Gaussian blur used');

%Gaussian Blurring Function

    sigblurx=sigall;
    sigblury=sigall;

    for i=Nx/2-4:Nx/2+5
        for j=Ny/2-4:Ny/2+5
            hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);
        end
    end

hh1=hh1/sum(sum(hh1));
hh1=fftshift(hh1);
Hcubic=conj(fft2(hh1));

Nx=Nx-2*extras;
Ny=Ny-2*extras;
    
hh=fftshift(hh);


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
q(Nx,1)=1;
q(1,Ny)=1;


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


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generate blurred data
%ssigma=sigma^2;
%n=(sigma)*randn(N/decFactor);
%g=imresize(real(ifft2(H.*fft2(f))),1.0/decFactor)+n;
name = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_LR_bicubic/X%d/0%dx%d.png',decFactor,imNumber,decFactor);

g = im2double(rgb2gray(imread(name)));

%g=g_(60/2:min([Nx Ny])/2+60/2-1,60/2:min([Nx Ny])/2+60/2-1);
%g=g_;
%%%%%%%%%%%%%%%%%%%%%%
%PSNRest=10*log10( 1/(ssigma*255) )

  
    
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
   

   
   
%==================Start Iterations==========================

%%
% x1=imresize( g,  decFactor  );
% 
% 
% %  hhz=zeros(N);
% %  v=fspecial('gaussian', 25,0.5);
% %  hhz(N/2-12+1:N/2+12+1,N/2-12+1:N/2+12+1)=v;
% %  HHZ=fft2(fftshift(hhz));
% 
% DFTShift=ones(Nx,Ny);
% 
% sigma_method=1;
% sigma=sqrt(ssigma);
% P=1;
% Nx_low=Nx/decFactor;
% Ny_low=Ny/decFactor;
% Lall=zeros(Nx_low,Ny_low,P);
% Lall(:,:,1)=g;
% dec = zeros(P,2); %This should be zero
% psigma = sigma*ones(P,1);
% dall   = zeros(P-1,2);
% 
% 
% %EM using stationary prior
% x=x1;
% 
% a=1;
% dfac=decFactor;
% 
% 
%   
% for iter=1:0
% 
%     DTemp(:,:,1) = DFTShift(:,:,1)/psigma(1);
%     Y1=fft2(  imupsampleblur(Lall(:,:,1), decFactor, dec(1,1:2),  conj(Hcubic))   ).*  DTemp(:,:,1)/psigma(1)   ;
% 
%     for k=2:P
%       [dummy, DFTShift(:,:,k)]=dft_shift(zeros(Nx,Ny),dall(k-1,:),N);
%       DTemp(:,:,k)=DFTShift(:,:,k)/(psigma(k));
%       Y1=Y1+fft2(imupsampleblur(Lall(:,:,k), dfac, dec(k,1:2),conj(Hcubic))).*conj(DTemp(:,:,k)/psigma(k));
%     end
% 
%     xprev=x1;
%     
%     [sum1,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I] = traceCom(ones(Nx,Ny),DTemp.*conj(Hcubic),abs(Q).^2,a,1,Nx,Ny,P,abs(Q).^2);
% 
%     X = invertY(Y1,Nx,Ny,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I);
%     x1 = real(ifft2( X ));
%     a = (NN-1)/(sum(sum(real(ifft2( Q.*X )).^2))+sum1)
% 
%    % One sigma for each image
%     if sigma_method==1
%      
%         s1 = real(HCfgCom(ones(Nx,Ny),DFTShift.*conj(Hcubic),A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I,Nx,Ny,P));
%         s2=0;
%  
%         for k=1:P
%             s2 = s2 + norm(Lall(:,:,k)-imdownsampleblur(x1,dfac,dec(k,1:2),conj(Hcubic).*DFTShift(:,:,k)),'fro')^2;
%         end
%         
%         ssigma=(s1+s2)/(P*NN/4)
%         psigma(:)=sqrt(ssigma);
%     end
%     
%     %Multiple sigma
%     if sigma_method==2
% 
% 
%         for k=1:P
% 
%             TEMP=ones(N,N,P);
%             for m=1:P
%                 TEMP(:,:,m)=conj(Hcubic).*DTemp(:,:,k)/sqrt(P);
%             end
%             
%             s1 = real(HCfgCom(ones(N,N),TEMP,A1I,A2I,B1I,B2I,C1I,C2I,D1I,D2I,N,P))
%             s2 = norm(  Lall(:,:,k)-imdownsampleblur(x1,dfac,dec(k,1:2),conj(Hcubic).*DFTShift(:,:,k)),'fro'   )^2
%             psigma(k)=sqrt((s1+s2)/(NN/4));
% 
%         end
%         psigma.^2
%     end
% 
% 
%     for k=2:P
%         [dall(k-1,:) DFTShift(:,:,k)] = shift_Newton(  fft2(Lall(:,:,k)   ),fft2(x1),dall(k-1,:),H(:,:,k),N);
%     end
% 
% end


g2=imresize(g,decFactor);


%ISNR_Bayes=10*log10( norm(f(20:Nx-20,20:Ny-20)  - (g2(20:Nx-20,20:Ny-20)),'fro')^2/norm( f(20:Nx-20,20:Ny-20) - x1(20:Nx-20,20:Ny-20),'fro')^2  )


name = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_HR/0%d.png',imNumber);

f_= im2double(rgb2gray(imread(name)));

fcolor=im2double(imread(name));


[Nx, Ny]=size(f_);

f=f_(1:Nx,1:Ny);



name = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_LR_bicubic/X%d/0%dx%d.png',decFactor,imNumber,decFactor);

g = im2double(rgb2gray(imread(name)));

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
        for j=1:Ny
            hh(i,j)=exp(-(i-Nx/2-1).^2/sigblurx).*exp(-(j-Ny/2-1).^2/sigblury);
        end
    end

%  hh=zeros(N);
%  v=fspecial('gaussian', 35,siggauss);
%  hh(N/2-12+1:N/2+12+1,N/2-12+1:N/2+12+1)=v;

  
else
    disp('Wrong blurtype argument');return;
end


hh=hh/sum(sum(hh));
if blurtype~=3
    hh=fftshift(hh);
end

H=fft2(hh);
expNum=100;
coord=zeros(expNum,2);


%optcoord=zeros(Nx,Ny);


W=zeros(Nx,Ny);
boundary=20;
W(boundary:Nx-boundary+1,boundary:Ny-boundary+1)=1;

indx2=0;
indx=1;
optcoord=zeros(Nx,Ny);
optcoord(Nx/2+1,Ny/2+1)=1;

for j1=1:2:9
    for k1=-j1:round(sqrt(j1)):j1+round(sqrt(j1))-1
        for l1=-j1:round(sqrt(j1)):j1+round(sqrt(j1))-1
        
  %            for k1=-j1:j1
  %              for l1=-j1:j1
             % if(k1~=0&&l1~=0)
            %(k1>0 | (k1==0 & l1>0))&
            
            if sqrt(k1^2+l1^2)<16  && (  k1>0  ||   (k1==0 && l1<0))  &&   (abs(k1)  ==   j1  ||   abs(l1)  ==   j1   )%(k1~=0 | l1~=0)%sqrt(k1^2+l1^2)<6&
            %   if (abs(k1)==j1 | abs(l1)==j1)&(k1~=0 | l1~=0)
                % l=l1
                %dirDiff(indx)  =   (norm( (W.*abs( fcolor(:,:,1)-circshift(fcolor(:,:,1), [k1 l1]  )   )).^2, 'fro' )  +   norm( (W.*abs( fcolor(:,:,2)-circshift(fcolor(:,:,2), [k1 l1]) )  ).^2, 'fro' )  +   norm( (W.*abs( fcolor(:,:,3)-circshift(fcolor(:,:,3), [k1 l1])  ) ).^2, 'fro' ));


                sigma1 = 1;
                sigma2 = 1;
                scale1 = 1;
                scale2 = 1;
                sigma1 = scale1*sigma1;
                sigma2 = scale2*sigma2;
                coord(indx,1)=mod(k1,Nx);
                coord(indx,2)=mod(l1,Ny);
                
                Theta = atan(l1/k1);
 
                hh=zeros(Nx,Ny);

                sigall=1;
                %Gaussian Blurring Function
                sigblurx=sigall*0.17;
                sigblury=sigall*0.17;
                for i=1:Nx
                    for j=1:Ny
                        hh(i,j)=exp(-(abs(i-Nx/2-1).^3.0)*sigblurx).*exp(-(abs(j-Ny/2-1).^3.0)*sigblury);
                    end
                end

                hh=imrotate(hh,Theta*180/(pi),'crop');
                %figure,imagesc(hh);

                S{indx}=1;
                hh=S{indx}*(hh)/sum(sum(hh));
                HQ{indx}=fft2(fftshift(hh));
  
                %if k1==0 && l1==0
                optcoord(l1+Nx/2+1,k1+Ny/2+1)=-1;

                %end
                    
                    
            indx=indx+1;
            end

        end
    end
end




expNum=indx-1;




for exper=1:expNum
    nu=0.000001;
    c2{exper}=1;
    rw{2*exper-1}(1)=coord(exper,1);
    rw{2*exper-1}(2)=coord(exper,2);
    Z{exper}=ones(Nx,Ny)/expNum;
    E{exper}=zeros(Nx,Ny)/expNum;
    A{exper}=zeros(Nx,Ny)/expNum;
    pof{exper}=1/expNum;
    Qp=zeros(Nx,Ny);
end;


g1__  =   im2double(imread(name));

g__  =   g1__(1:Nx/decFactor,1:Ny/decFactor,:);

g_r=g__(:,:,1);
g_g=g__(:,:,2);
g_b=g__(:,:,3);

x_r=zeros(Nx,Ny,1);
x_g=zeros(Nx,Ny,2);
x_b=zeros(Nx,Ny,3);

g2=imresize(g__,decFactor);
%g2_=imsharpen(g2,'Radius',3,'Amount' ,3);
x_r=g2(1:Nx,1:Ny,1);
x_g=g2(1:Nx,1:Ny,2);
x_b=g2(1:Nx,1:Ny,3);

g1=zeros(Nx,Ny);





ssigma=0.000001;

for iter=1:5

     ZALL=zeros(Nx,Ny);

     for exper=1:expNum

        J=(S{1}/2+nu/2);

        qa=zeros(Nx,Ny);
        qa(1,1)=1;
        qa(coord(exper,1)+1,coord(exper,2)+1)=-1;

        q1=real(ifft2(fft2(qa))).^2;
        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
      
        E{exper}   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2  + real(ifft2(fft2(q1).*(fft2(Qp))))) .* conj(HQ{exper})  ));
        % E{exper}   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* HQ{1}  ));

        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
        A{exper}  =   (nu+S{1})./(c2{exper}*E{exper}+nu)   ;
        Z{exper}  =   exp(  J*log(A{exper})-A{exper}.*(c2{exper}*E{exper}/2+nu/2)   );


        %Z{exper}=pof{exper}.*(nu+c2{exper}*E).^(-J);
        %B{exper}=A{exper};
        ZALL=ZALL+Z{exper};

    end

    for exper=1:expNum
        Z{exper}=Z{exper}./ZALL;
        %Z{exper}(find(Z{exper}<0.05))=0.0;
        B{exper}=A{exper}.*Z{exper};
        B{exper}=  c2{exper}  *   real(ifft2(  fft2(  B{exper}   )  .*  ( HQ{exper})   ))  /   S{1};
        rw{2*exper}  =   gpuArray(  ssigma*B{exper}   );
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
%         %c1{exper}=bisection(0.00000812,100,B{exper}./ZALL,Z{exper},c1{exper},S,N);
%     end

    Ainv=single( real(ifft2(zeros(Nx,Ny))) );
%    xprev=x;
    
    Qp=sum(sum( real(ifft2( H )).*real(ifft2( conj(H) )) ))/ssigma;

    for exper=1:expNum

        %pof{exper}=exp( psi(Z{exper}+10000) -psi(A0) );%(Z{exper}+0.0812)./ A0;
        qa=zeros(Nx,Ny);
        qa(1,1)=1;

        qa(coord(exper,1)+1,coord(exper,2)+1)=-1;
        q1=real(ifft2(fft2(qa))).^2;
        Qp=Qp+real(ifft2(fft2(q1).*(fft2(B{exper}))));
    end

    Qp=1./Qp;

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

    % x=reshape(cgLanczos(@Amat,zeros(NN,1), b(:), iw, rw, 1, 0, 550, 0.0812),N1,N2);

    %toc
%     iter
%     if iter == 3 || iter==6
%         x=x+0.01*randn(N);
%     end
    itnum=0;
%         for exper=1:expNum
%             Qp{exper}=zeros(N);
%         end



 

        for iter1=1:itnum
            randn('state',iter1);
            d=randn(Nx,Ny);
            
            %e  = cudaSymLanczosBlockMatch(single(zeros(N)),single(real(ifft2(H))/sqrt(ssigma)),expNum,coord',B,single(ones(N)),single(d),ssigma,N,N,R,0.0812,150) ;
            iw=abs(H).^2/ssigma;
            b=d;
            tic;e=reshape(cgLanczos(@Amat,zeros(NN,1), b(:), iw, rw, 1, 0, 1050, 0.000812),N,N);toc

            for exper=1:expNum
                Qp{exper} = Qp{exper} +  (e-circshift(e,rw{2*exper-1})) .*  (d-circshift(d,rw{2*exper-1}))/itnum;
            end


        end




    if itnum>0

        for iter1=1:itnum
            randn('state',iter1);
            d=randn(Nx,Ny);
            
            %e  = cudaSymLanczosBlockMatch(single(zeros(N)),single(real(ifft2(H))/sqrt(ssigma)),expNum,coord',B,single(ones(N)),single(d),ssigma,N,N,R,0.0812,150) ;
            iw=abs(H).^2/ssigma;
            b=d;
            tic;e=reshape(cgLanczos(@Amat,zeros(NN,1), b(:), iw, rw, 1, 0, 1050, 0.000812),N,N);toc

            for exper=1:expNum
                Qp{exper} = Qp{exper} +  (e-circshift(e,rw{2*exper-1})) .*  (d-circshift(d,rw{2*exper-1}))/itnum;
            end

        end
    end


    iw{1}=gpuArray(H);
    iw{2}=decFactor;
    iw{3}=Hcubic;
    %imresize(  g,   decFactor);
    %g1(2:decFactor:Nx, 2:decFactor:Ny)=g;

    %g1=zeros(Nx,Ny);
    %g1(1:decFactor:Nx,1:decFactor:Ny)=g_r;
    %g1 = real(ifft2( fft2(g1) .* Hcubic));%imresize(g_r,decFactor);
    g1 = myimresizeTranspose(g_r,Nx,Ny,Hcubic,decFactor,4) ;
  
    x_rprev=x_r;
    b = gpuArray(g1);
    tic;   x_r=reshape( cgLanczos( @Amat,gpuArray(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 150, 0.0000001), Nx, Ny ); toc

   % g1 = imresize(g_g,decFactor);
    %g1=zeros(Nx,Ny);
    %g1(1:decFactor:Nx,1:decFactor:Ny)=g_g;
    %g1 = real(ifft2( fft2(g1) .* Hcubic));%imresize(g_r,decFactor);
    g1 = myimresizeTranspose(g_g,Nx,Ny,Hcubic,decFactor,4) ;

    b = gpuArray(g1) ;

    x_gprev=x_g;

    tic; x_g=reshape(cgLanczos(@Amat,gpuArray(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 150, 0.0000001), Nx, Ny); toc

    %g1 = imresize(g_b,decFactor);

    %g1=zeros(Nx,Ny);
    %g1(1:decFactor:Nx,1:decFactor:Ny)=g_b;
    %g1 = real(ifft2( fft2(g1) .* Hcubic));%imresize(g_r,decFactor);
    g1 = myimresizeTranspose(g_b,Nx,Ny,Hcubic,decFactor,4) ;
    
    b = gpuArray(g1);

    x_bprev=x_b;
    tic;x_b=reshape(cgLanczos(@Amat,gpuArray(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 150, 0.0000001), Nx, Ny); toc

    %figure,imagesc(x_b)
    %ISNR_Bayes=10*log10(norm(f-imresize(g,decFactor),'fro')^2/norm(f-x,'fro')^2)


    %g__2=imresize(g,decFactor);
    %ISNR_Bayes=10*log10(norm(fcolor(20:Nx-20,20:Ny-20)-g__2(20:Nx-20,20:Ny-20),'fro')^2/(norm(fcolor(20:Nx-20,20:Ny-20,1)-x_r(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,2)-x_g(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,3)-x_b(20:Nx-20,20:Ny-20),'fro')^2))

    boundary=decFactor+7;
    
    
%     x_r(1:boundary,Ny-boundary:Ny)=g2(1:boundary,Ny-boundary:Ny,1);
%     x_g(1:boundary,Ny-boundary:Ny)=g2(1:boundary,Ny-boundary:Ny,2);
%     x_b(1:boundary,Ny-boundary:Ny)=g2(1:boundary,Ny-boundary:Ny,3);
%     
%     x_r(1:boundary,1:boundary)=g2(1:boundary,1:boundary,1);
%     x_g(1:boundary,1:boundary)=g2(1:boundary,1:boundary,2);
%     x_b(1:boundary,1:boundary)=g2(1:boundary,1:boundary,3);
%     
%     x_r(Nx-boundary:Nx,1:boundary)=g2(Nx-boundary:Nx,1:boundary,1);
%     x_g(Nx-boundary:Nx,1:boundary)=g2(Nx-boundary:Nx,1:boundary,2);
%     x_b(Nx-boundary:Nx,1:boundary)=g2(Nx-boundary:Nx,1:boundary,3);
%     
% 
%     x_r(Nx-boundary:Nx,Ny-boundary:Ny)=g2(Nx-boundary:Nx,Ny-boundary:Ny,1);
%     x_g(Nx-boundary:Nx,Ny-boundary:Ny)=g2(Nx-boundary:Nx,Ny-boundary:Ny,2);
%     x_b(Nx-boundary:Nx,Ny-boundary:Ny)=g2(Nx-boundary:Nx,Ny-boundary:Ny,3);

   ISNR=10*log10( ( norm( fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1)  -  g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),'fro')^2   +norm( fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,2)  -  g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1,2),'fro')^2+norm( fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,3)  -  g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1,3),'fro')^2)/(  norm(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1)-x_r(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2+norm(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,2)-x_g(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2+norm(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,3)-x_b(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

    PSNR=10*log10( (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)*3/(   norm(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1)-x_r(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro'    )^2  +   norm(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,2)  -   x_g(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2+norm(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,3)-x_b(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))


   %ISNR_Bayes=10*log10(norm(f(20:Nx-20,20:Ny-20)  -  g_(20:Nx-20,20:Ny-20),'fro')^2/norm(f(20:Nx-20,20:Ny-20)-x(20:Nx-20,20:Ny-20),'fro')^2)

   SSIMres=(ssim(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x_r(boundary:Nx-boundary+1,boundary:Ny-boundary+1)))+ssim(fcolor(boundary:Nx-boundary+1,boundary:Ny-boundary+1,2),gather(x_g(boundary:Nx-boundary+1,boundary:Ny-boundary+1)))+ssim(fcolor(boundary:Nx-boundary,boundary:Ny-boundary+1,3),gather(x_b(boundary:Nx-boundary,boundary:Ny-boundary+1))))/3.0;


%           
%         for exper=1:expNum
%             
%             E1 = ( ( x_r-circshift(x_r,rw{2*exper-1}) ).^2  + ( x_g-circshift(x_g,rw{2*exper-1}) ).^2 +( x_b-circshift(x_b,rw{2*exper-1}) ).^2  );
%             %A1 = (nu+1)./(real(ifft2( fft2( E1 ).*conj(HQ{exper}))) +nu);
%             %Z{exper}=Z{exper}./ZALL;
%             %              %Z{exper}=(gamma(J)/gamma(c1{exper}/2))*(c2{exper}/c1{exper})^(S/2)*(1+c2{exper}*E/c1{exper}).^(-J);
%             %              %Z(find(real(ifft2(fft2(Z.*HQ{exper})))<0.001))=0.000000001;
%             %              %B{exper}=Z{exper}.*A;
%             U = real(ifft2( conj(fft2( (fftshift( (A{exper}.*Z{exper}))))) .* fft2(E1) )) ;
% 
%             %R=5;
% 
%             %h2=exp( - real(ifft2( conj(fft2(  A1.*Z{exper})) .* fft2(fftshift(E1)) )) );
%             %res  =   real(ifft2( conj(fft2(A{exper}.*Z{exper})) .* fft2(fftshift(E1))));
%             %rmax=max(max(res));
%             %h2=exp(-real(ifft2( conj(fft2(A{exper}.*Z{exper})) .* fft2(fftshift(E1)))) +rmax );
%             Nxy=size(find(Z{exper}>0.1));
%             h2=exp(-0.5*U/Nxy(1));
% 
% %            m=zeros(Nx,Ny);
% %            for i=Nx/2-R:Nx/2+R+1
% %                for j=Ny/2-R:Ny/2+R-1
% %                    U=0.5*circshift( E1 ,[i-Nx/2, j-Ny/2]).*Z{exper};
% %                    %m(i,j)  =   gather( sum(sum(  (  0.5*(    exp(-U)    )      )     )          )        );
% %                    m(i,j)  =   gather( sum(sum(  (  0.5*(    U.^(-(Z{exper}/2+1/2)    )      )     )          )        ));
% %                    %figure(1),imagesc(h),colormap('gray');
% %                end
% %            end
% 
% 
%             %m=m/sum(sum(m));
%             h2=fftshift(h2)/sum(sum(h2));
%             HQ{exper}=(fft2(h2));
% 
%         end

end
 

%save(sprintf('file%dwithISNR%d',myinput1),'ISNR_Bayes');
