function [PSNR, ISNR, SSIMres]=SRviaDenoisingYCbCx3Gauss_26_1_2020(pathHR,pathLR,decFactor)


%Changes: for 3x decFactors


%clear all;
%myinput1=1;
%imNumber=2;
%randn('seed',0);


decFactor=3;

%pathHR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_train_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);
%pathLR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);
%pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);

%pathHR = sprintf('/home/gchantas/Downloads/Set14x1_%d/image_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Downloads/Set14x%d/image_%03d.png',decFactor,imNumber);

%pathHR = sprintf('/home/gchantas/Downloads/Set5x1_%d/image_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Downloads/Set5x%d/image_%03d.png',decFactor,imNumber);



%name ='barbara.png';
%name ='Lena512.png';
%name ='boat.png';
%name =
%sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
extras=0;%expand by $extras pixels the boundary of the image by replicating the bounds


f_= im2double(imread(pathHR));

[Nx, Ny,dummy]=size(f_);

f=f_(1:Nx,1:Ny,1);

[Nx, Ny]=size(f);


padx=2*ceil(Nx/2)-Nx
pady=2*ceil(Ny/2)-Ny


f=padarray(f,[padx pady],'post');
[Nx, Ny]=size(f);

g = im2double(  imread(pathLR)   );
%g=g;
g2=imresize( g, 3);
g2=padarray(g2,[padx pady],'post');

  %Nx=min(1024,Nx);
  %Ny=min(1024,Ny);
%f=f_(100:Nx+99,100:Nx+99);
%f=f_(1:Nx,1:Ny);

%g(1:min(Nx,Nx1),1:min(Ny,Ny1))=f1(500:499+min(Nx,Nx1),1:min(Ny,Ny1),1);



% 
%     sigall=1;
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
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generate Bicubic interpolator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hh1=zeros(Nx,Ny);
sigall=1/2;
 
%Bicubic separable filter

%     sigblurx=sigall;
%     sigblury=sigall;
% 
%     for i=round(Nx/2)-12:round(Nx/2)+13
%         for j=round(Ny/2)-12:round(Ny/2)+13
%          %     hh1(i,j)=sigblurx*cubic(-(i-Nx/2+1/4)*sigblurx-1/16)*sigblury*cubic(-(j-Ny/2+1/4)*sigblury-1/16);
%                        
%            %For decFactor=2
%        hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);
% 
%         end
%     end
% 
% hh1=hh1/sum(sum(hh1));
% hh1=fftshift(hh1);
% Hcubic=(fft2(hh1));

    

%norm(myimresize(f,Nx,Ny,Hcubic,decFactor,extras)-imresize(f,1/decFactor),'fro')^2/(Nx*Ny)


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

Qp=abs(fft2(Q)).^2;



alpha=1;
%ssigma=10^(-8+decFactor/2); previous runs 1/1/2019
ssigma=10^(-7)/2;

x=g2;

[x,alpha,ssigma,Inv,HDDH]=StatSR(x,imresize(g2,0.5),1,ssigma,Nx,Ny,2,f);
alpha=alpha/4

%x(30:Nx-30,30:Ny-30)=(xprev(30:Nx-30,30:Ny-30));

'Stationary'
norm(f-x,'fro')^2/(Nx*Ny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f_= im2double(imread(pathHR));

[Nx, Ny,dummy]=size(f_);

f=f_(1:Nx,1:Ny,1);

[Nx, Ny]=size(f);

g = im2double(  imread(pathLR)   );
%g=g;
g2=imresize(g, decFactor);
hNx=floor(Nx/2);
hNy=floor(Ny/2);
x=x(1:Nx,1:Ny);
% 
% 
% psInv=abs(Hcubic).^2;
% psInv=psInv(1:Nx/2,1:Ny/2)+psInv(Nx/2+1:Nx,Ny/2+1:Ny)+psInv(1:Nx/2,Ny/2+1:Ny)+psInv(Nx/2+1:Nx,1:Ny/2);
% 
% indices=find(Hcubic<0.000000001);
% notindices=find(psInv>=0.00001);
% psInv(psInv<100*eps)=0;
% 
% psInv(psInv~=0)=1./psInv(psInv~=0);
% %psInv(notindices)=1./psInv(notindices);
% %psInv(indices)=1;
% Hcubic(indices)=1;
% x1=real(ifft2(fft2(g).*(psInv)));
% x1_=zeros(Nx,Ny);
% x1_(1:decFactor:Nx,1:decFactor:Ny)=x1;
% x2=decFactor^2*real(ifft2(conj(Hcubic).*fft2(x1_)));
% norm(f-x2,'fro')^2/(Nx*Ny)
% 
% 
% 'Pseudo'
% norm(f-x2,'fro')^2/(Nx*Ny)
% return
%%
Nx=Nx+2*extras;
Ny=Ny+2*extras;

hh1=zeros(Nx,Ny);
sigall=1/3;
 


%Bicubic separable filter

    sigblurx=sigall;
    sigblury=sigall;

    for i=floor(Nx/2)-12:floor(Nx/2)+13

        for j=floor(Ny/2)-12:floor(Ny/2)+13
           %For decFactor=4
            % hh1(i,j)=sigblurx*cubic(-(i-hNx+1/4)*sigblurx-1/16)*sigblury*cubic(-(j-hNy+1/4)*sigblury-1/16);
           %For decFactor=2
           %hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);
            %hh1(i,j)=sigblurx*cubic(-(i-floor(Nx/2)+1/3)*sigblurx)*sigblury*cubic(-(j-floor(Ny/2)-1+1/3)*sigblury);
            %hh1(i,j)=sigblurx*cubic(-(i-floor(Nx/2)-1+1/2)*sigblurx*3/2)*sigblury*cubic(-(j-floor(Ny/2)-1+1/2)*sigblury*3/2);
          
          hh1(i,j)=sigblurx*cubic(-(i-round(Nx/2))*sigblurx)*sigblury*cubic(-(j-round(Ny/2))*sigblury);
          %For decFactor=3
          %  hh1(i,j)=sigblurx*cubic(-(i-floor(Nx/2)+1/3)*sigblurx-1/9)*sigblury*cubic(-(j-floor(Ny/2)+1/3)*sigblury-1/9);
            
            
        end

    end

hh1=hh1/sum(sum(  hh1   ));
hh1=fftshift(hh1);
Hcubic=fft2(hh1);

Nx=Nx-2*extras;
Ny=Ny-2*extras;
    

%Qp=zeros(N);
   
%ISNR_Bayes=10*log10( norm(f(20:Nx-20,20:Ny-20)  - (g2(20:Nx-20,20:Ny-20)),'fro')^2/norm( f(20:Nx-20,20:Ny-20) - x1(20:Nx-20,20:Ny-20),'fro')^2  )

norm(myimresize(f,Nx,Ny,Hcubic,decFactor,extras)-imresize(f,1/decFactor),'fro')^2/(Nx*Ny)



expNum=400;
coord=zeros(expNum,2);

indx2=0;
indx=1;
optcoord=zeros(Nx,Ny);
optcoord(hNx+1,hNy+1)=1;

            for exper=1:expNum
                nu{exper}=4.0001;
                c2{exper}=5000;
                Z{exper}=ones(Nx,Ny)/expNum;
                A{exper}=zeros(Nx,Ny)/expNum;
                pof{exper}=1/expNum;
                Qp=zeros(Nx,Ny);
                E{exper}=zeros(Nx,Ny)/expNum;
            end;


            %x=g2(1:Nx,1:Ny);
            



indx2=0;
indx=1;
optcoord=zeros(Nx,Ny);
optcoord(hNx+1,hNy+1)=1;
%x=g2;

%corrwin=15;
% Fmask=zeros(Nx,Ny);
% Fmask(1:10,1:2)=1;
% Fmask(1:2,1:10)=1;
% 
% Fmask2=zeros(Nx,Ny);
% Fmask2(1:10,1:2)=1;
% Fmask2(1:10,Ny)=1;
% Fmask2(1:2,1:10)=1;
% Fmask2(Nx,1:10)=1;
% Fmask2(1,1)=0;
% fcorr=fftshift(real(ifft2(abs(fft2((g2)).^2))));

% 
% [R,angles]=radon(fcorr);
% [x1,y1]=sort(R(:),'descend');
% [k1, l1]= ind2sub(size(R),y1(1));
% 
% Fmask=fftshift(Fmask)+imrotate(fftshift(Fmask2),45,'crop');
% Fmask=imrotate(Fmask,180-l1,'crop');
% %figure,imagesc((Fmask));
% 
Fmask=zeros(Nx,Ny);
Fmask(hNx+1:2:hNx+1+90,hNy-46:2:hNy+46)=1;
Fmask=Fmask+circshift(Fmask,[1 1]);%Fmask(Nx/2+1:2:Nx/2+1+50,Ny/2-25:1:Ny/2+25)=1;


Fmask(hNx+1:hNx+8,hNy-4:hNy+4)=1;


Fmask(hNx+1,hNy+1:Ny)=0;

g3=zeros(Nx*2,Ny*2);
g3(1:Nx,1:Ny)=g2;
%g3(1:Nx,1:Ny)=g3(1:Nx,1:Ny);
%fcorr=fftshift(real(ifft2(abs(fft2((g2)).^2))));
fcorr=fftshift(real(ifft2(abs(fft2((g3-mean(mean(g3)))).^2))));
fcorr=fcorr(floor(hNx)+1:floor(Nx/2*3),floor(Ny/2)+1:floor(Ny/2*3));
fcorrMax=fcorr(hNx+1,hNy+1);

fcorr=fcorr.*(Fmask);


%maxExp=sum(sum(Fmask));
%fcentral=fcorr(Nx/2-corrwin:Nx/2+corrwin+1,Ny/2-corrwin:Ny/2+corrwin+1);
%fcorrsorted=sort(fcentral(:));

tic;[fcorrsorted,I]=sort(fcorr(:),'descend');toc
[k1, l1]= ind2sub([Nx,Ny],I(:)); 
%fcorrsorted(1:100,ind2sub([Nx,Ny],1:100));
%optcoord(k1(1:100), l1(1:100))=-1;
%[k1, l1]=find(fcorrsorted(Nx-10:Nx,Ny-10:Ny));
%[k1, l1]=find(fcorr(1:1:Nx,1:1:Ny)>median((fcentral(:))));
% RFCorr=radon(fcorr);
% [Nr1, Nr2]=size(RFCorr);
% [maxangle, MI] =max(RFCorr(:));
% [k1, l1]= ind2sub([Nx,180],I(:));




expNum=1;
indx=0;

for j1=1:1:280

    indx=j1;

    if mod(k1(indx)-hNx-1,Nx)~=0 || mod(l1(indx)-hNy-1,Ny)~=0
    
        coord(expNum,1)=mod(k1(indx)-hNx-1,Nx);
   
        coord(expNum,2)=mod(l1(indx)-hNy-1,Ny);

        optcoord(k1(indx), l1(indx))=-1;

        expNum=expNum+1;
    else
       [ mod(k1(indx)-hNx-1,Nx), mod(l1(indx)-hNy-1,Ny)]
    end
    
    if  fcorr(k1(indx),l1(indx))<fcorrMax*0.60
     break;
    end

end

expNum=expNum-1



            %Kernel for patch congruency modeling
%             sigblurx=sigall*0.5;
%             sigblury=sigall*0.5;
%              for i=1:Nx
%                 for j=1:Ny
%                     hh(i,j)=exp(-(abs(i-floor(Nx/2)-1).^3.0)*sigblurx).*exp(-(abs(j-floor(Ny/2)-1).^3.0)*sigblury);
%                 end
%              end
 %Gaussian Blurring Function
            sigblurx=0.25;
            sigblury=0.25;
             for i=1:Nx
                for j=1:Ny
                    hh(i,j)=exp(-(abs(i-hNx-1).^3.0)*sigblurx).*exp(-(abs(j-hNy-1).^3.0)*sigblury);
                end
             end

            S{1}=1;
            hh=S{1}*(hh)/sum(sum(hh));
            HQ{1}=fft2(fftshift(hh));

            %x=g2(1:Nx,1:Ny);
            
% 
% g3=zeros(Nx,Ny);
% g3(1:Nx,1:Ny)=f;
% 
% 
% for iter=1:30
% 
% 
% hh=fftshift((real(ifft2(abs(fft2((g3-mean(mean(g3)))).^2)))))/(0.5*sum(sum(g3-mean(mean(g3))).^2  )  );
% hh(hh<0)=0;
% 
% %hh(1:Nx,1:Ny)=hh_(Nx/2+1:Nx*3/2,Ny/2+1:Ny*3/2);
% 
% Imask2=zeros(Nx,Ny);
% Imask2(Nx/2+1-30:Nx/2+1+30,Ny/2+1-30:Ny/2+1+30)=1;
% hh=Imask2.*hh;
% 
% 
% hh(hh~=0)=exp(-1./(hh(hh~=0)));
% S{1}=1;
% hh=S{1}*(hh)/sum(sum(hh));
% HQ{1}=fft2(fftshift(hh));
% g3=real(ifft2( fft2(f).*HQ{1}));
% end


% 
% % 
% % 
%   for exper=1:expNum
%       E{exper}  =  gpuArray( real(  ifft2( fft2(  (f-circshift(f, coord(exper,:)  )).^2   )  .* conj(HQ{1})   )));
%       %  E=real(  ifft2( fft2( ( g2-circshift(g2, coord(exper,:))  )  .* conj(HQ{1})  )));
% %       
% %       %  VE{exper}   =  real(  ifft2( fft2( ( g2-circshift(2, coord(exper,:)) -E ).^2 ) .* conj(HQ{1})  ));
%   end
% %  
% %  
% %     
%   [loglikelihood,nu,c2]=trainPCMixtureOneL(expNum,c2,nu,E,S);
%  

%[loglikelihood,nu,c2]=trainPCMixture(expNum,c2,nu,coord,E,S);
 
% [nu,c2]=learnPCMixture(expNum,x,c2,nu,Hcubic,HQ{1},ssigma,S,rw,decFactor,extras)
% 
 

% % %    
% % 
% % 
% %       disp('Nus')
% %       return
%  
%       nu{1:expNum}
         
%       coordbackup=coord;
% 
% coord=coordbackup;
% coord2=zeros(expNum,2);
% indx=1;
% 
% for exper=1:expNum, if nu{exper}/(c2{exper}.*(nu{exper}-2))<0.00005,  coord2(indx,:)=coord(exper,:); nu{indx}=nu{exper};c2{indx}=c2{exper}; indx=indx+1;end; end;
% coord=coord2;
% 
% 
% 
% expNum=indx-1;
 
 %return;
clear E;


        for exper=1:expNum
            %nu=0.00001;
            %c2{exper}=1.01;
            rw2{exper}(1)=coord(exper,1);
            rw2{exper}(2)=coord(exper,2);
            Z{exper}=ones(Nx,Ny)/expNum;
            A{exper}=zeros(Nx,Ny)/expNum;
            %pof{exper}=1/expNum;
            Qp=zeros(Nx,Ny);
            PW{exper}=ones(Nx,Ny)/expNum;
        end;
%x=g2(1:Nx,1:Ny);

        ssigma=10^(-7)/(1);%/log(alpha+1);
    
        
     %numbers indicating by how many pixels calculation of differences goes beyond the
     %image boundary
     
     for exper=1:expNum
         maxX(exper)=max(abs( mod((coord(exper,1)+hNx),Nx)-hNx))+1;
         maxY(exper)=max(abs( mod((coord(exper,2)+hNy),Ny)-hNy))+1;
     end
     
    maxmaxX=max(maxX);
    maxmaxY=max(maxY);
     
     

z=gpuArray2(g2);

for iter=1:200
tic

     ZALL=zeros(Nx,Ny);

     for exper=1:expNum

        J=(S{1}/2+nu{exper}/2);

        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;

         E   =  real(  ifft2( fft2( ( x-circshift(x, rw2{exper})  ).^2  + Qp+circshift(Qp, rw2{exper})) .* conj(HQ{1})  ));
      %  E   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* conj(HQ{exper})  ));

        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
        A  = ones(Nx,Ny);%  (nu{exper}+S{1})./(c2{exper}*E+nu{exper})   ;
        %Z{exper}  =   exp(  J*log(A{exper})-A{exper}.*(c2{exper}*E/2+nu{exper}/2)   );
        %Z{exper}  =   exp(  log(c2{exper})/2+(S{1}/2+nu{exper}/2-1)*(log(A{exper})-log(J)+psi(J))-A{exper}.*(c2{exper}*E/2+nu{exper}/2) -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2)  );


       % Z=   ((gamma((nu{exper}+S{1})/2)/gamma((nu{exper})/2)) *(c2{exper}/nu{exper})^0.5*(1+c2{exper}*E/nu{exper} ).^(-(nu{exper}+S{1})/2 ));%
        Z=exp(-0.5*c2{exper}*E);

        %c2prev{exper}=c2{exper};

        %if iter >1 
            %c2{exper}=Nx*Ny/sum(sum( (B{exper}/c2prev{exper}).*( x-circshift(x, rw{2*exper-1})  ).^2));
        %end

        %Z{exper}=(nu+c2{exper}*E).^(-J);
        rw{exper}  =   A.*Z;
        ZALL=ZALL+Z;

     end

   % for exper=1:expNum
    %   Z{exper}=Z{exper}./ZALL;
    %   PW{exper}=(Z{exper}+.1)/(1+.1);
       % pof{exper}=sum(sum(Z{exper}));
    %end
% Imask3=zeros(Nx,Ny);
% bmaxX=Nx/4+10;
% bmaxY=Ny/4+10;
%Imask3(Nx/2-bmaxX:Nx/2+1+bmaxX,Ny/2-bmaxY:Ny/2+1+bmaxY)=1;
    Qp=zeros(Nx,Ny);

    
    for exper=1:expNum
        %Z{exper}=Z{exper}./ZALL;
        %Z{exper}(find(Z{exper}<0.05))=0.0;
        B=rw{exper} ;
        B=B./ZALL;
        

      %   E   =  real(  ifft2( fft2( ( x-circshift(x, rw{2*exper-1})  ).^2  + Qp+circshift(Qp, rw{2*exper-1})) .* conj(HQ{1})  ));
       % c2prev=sum(sum(Z{exper}))*nu{exper}/sum(sum(B{exper}.*E));
       B=  c2{exper}  *   real(ifft2(  fft2(  B   )  .*   ( HQ{1})   ))  /   S{1};

      %B{exper}=Imask3.*fftshift(ifft2(B{exper}));
     % B{exper}=real(fft2(fftshift(B{exper})));
      % B{exper}(B{exper}<100)=0;
        %[ZeroIndX(exper) ZeroIndY(exper)]=find(B{exper}<.1);
        
        rw{exper}  =   gpuArray2(  B  );

      %  c2{exper}=c2prev;
      % disp(sprintf('c2(%d) %f', exper, c2{exper}));
       %Qp=Qp+B+circshift(B,rw{2*exper-1});%real(ifft2(conj(fft2(q1)).*(fft2(B{exper}))));
               
    end
        %disp 'time to B, E and Z'
%toc
    
%x=g2;

alpha=2000*sqrt(iter);%.^1.1;

ssigma=10^(-7);

x_prev=x;
    




for iter2=1:100

%     iw{1}=alpha;
%     iw{2}=decFactor;
%     iw{3}=ones(Nx,Ny);
%     iw{4}=extras;
    %imresize(  g,   decFactor);
    %g1(2:decFactor:Nx, 2:decFactor:Ny)=g;

    %g1=zeros(Nx,Ny);
    %g1(1:decFactor:Nx,1:decFactor:Ny)=g_r;
    %g1 = real(ifft2( fft2(g1) .* Hcubic));%imresize(g_r,decFactor);
    %g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;

%     b = gpuArray2(x)*alpha;
%     tic;   [ z, istop, itn, Anorm, Acond, rnorm, xnorm, D_r ]=cgLanczos( @DenoiseAmat,gpuArray(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 400, 10^(-10)); toc
% 
%  z=reshape(z,  Nx, Ny );

    % PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - z(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )
    %          if norm(z-z1,'fro')<0.0001
    %              'Number of iteration vardenoise stopped'
    %              iter3
    %              break
    %          end
    

    z1=gpuArray2(z);

    
    for exper=1:expNum
         maxX(exper)=max(abs( mod((coord(exper,1)+Nx/2),Nx)-Nx/2))+1;
         maxY(exper)=max(abs( mod((coord(exper,2)+Ny/2),Ny)-Ny/2))+1;
     end
     
    maxX=max(maxX);
    maxY=max(maxY);
    %maxX=30;
    %maxY=30;

    for exper=1:expNum
        rw{exper}=padarray(rw{exper},[maxX maxY],'circular');
    end

    %tic

    for iter3=1:1

        zprev=z;

        [z,Ball]=VariationalDenoise(padarray(z1,[maxX maxY],'circular'),padarray(x,[maxX maxY],'circular'),rw,Nx,Ny,coord,expNum,alpha,maxX,maxY);
            norm(z-z1,'fro')

        if norm(z-z1,'fro')^2/(Nx*Ny)<10^(-7)
            norm(z-z1,'fro')
            disp('z evaluation converged')
            break;
        end;

        z1=z;

    end



    for exper=1:expNum
        rw{exper}=rw{exper}(maxX+1:Nx+maxX,maxY+1:Ny+maxY);
    end
 
   % toc

    'time to var x'
   
    Qp=1./(alpha+Ball);
    
   
   %PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - z(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )

   
   %z=reshape(z,  Nx, Ny );
   %g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
 
   %b = gpuArray2(g1+alpha*ssigma*z);
   
   % x1=real(ifft2( HF.* padarray(b,[extras extras],'circular')));
   % x=x1(extras+1:Nx-extras,extras+1:Ny-extras);

  iw{1}=alpha*ssigma;
  iw{2}=decFactor;
  iw{3}=Hcubic;
  iw{4}=extras;

  g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras);

  b = gpuArray2(g1+alpha*ssigma*z);
  tic
  xprev=x;
  x  =   cgLanczos(  @EasySRAmat,gpuArray2(zeros(Nx,Ny))   , b(:), iw, rw, 1, 1, 1000, 0.00000001); 
  toc
  
  x=reshape(x,  Nx, Ny );

%alpha=1;
%z_=padarray(z,[extras extras],'symmetric');
%Z_=fft2(z_);
%g3=imresize(g,2);

%tic
% 
% decFactor=decFactor^2/2;
% 
% Ia=makeDiagMatCells(  alpha*ssigma*ones(Nx+2*extras,Ny+2*extras),Nx+2*extras,Ny+2*extras,decFactor*2,0   );%%diagonal identity matrix (no blocks)
% C2=addBlockMat(C,Ia,decFactor*2);
% Inv=invertBlockMatrix(C2,decFactor*2);
% 
% G3=makeBlockVec( G_ ,Nx+2*extras,Ny+2*extras,decFactor,0);%applyBlockDiag(HcubBlockConj,G2_,decFactor*2);
% bVec=makeBlockVec(  ssigma*alpha*Z_,Nx+2*extras,Ny+2*extras,decFactor,0);
% bVec2=addBlockVec(G3,bVec,decFactor*2);
% 
% xblock=applyBlockDiag(Inv,bVec2,decFactor*2);
% 
% vec=deblockVec(xblock,Nx+2*extras,Ny+2*extras,decFactor*2,0);
% xprev=x;
% x=real(ifft2(vec));
% x=x(extras+1:Nx+extras,extras+1:Ny+extras);
% 
% decFactor=dcf;
%diagHHa=tageDiagBlockMat(Inv,Nx,Ny,decFactor*2,0);
%toc


   %D_r=reshape(D_r,  Nx, Ny );
   %ISNR_Bayes=10*log10(norm(fcolor(20:Nx-20,20:Ny-20)-g__2(20:Nx-20,20:Ny-20),'fro')^2/(norm(fcolor(20:Nx-20,20:Ny-20,1)-x_r(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,2)-x_g(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,3)-x_b(20:Nx-20,20:Ny-20),'fro')^2))

   boundary=decFactor;

   ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

   
   %g2 is the bicubic SR image
   %PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )

   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )

   SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));
   
% 
%   %  Nx=Nx+2*extras;
%   %  Ny=Ny+2*extras;
%     %Hcubic=fftshift(Hcubic);
%     Imask=ones(Nx,Ny);
%     Imask(1:decFactor:Nx,1:decFactor:Ny)=0;
%     Imask= not(Imask);
%     hcubatrous=fftshift(real(ifft2( Hcubic )));
%     hcubatrous=hcubatrous(extras+1:Nx+extras,extras+1:Ny+extras);
%     hcubatrous=hcubatrous.*Imask;
   
%  HF= hcubatrous.^2;%abs(Hcubic(1:Nx/2,1:Ny/2)).^2+abs(Hcubic(Nx/2+1:Nx,1:Ny/2)).^2+abs(Hcubic(1:Nx/2,Ny/2+1:Ny)).^2+abs(Hcubic(Nx/2+1:Nx,Ny/2+1:Ny)).^2;
%  
%     %'alpha denumerator'
%     TrH= sum(sum( HF./(HF+ssigma*alpha*ones(size(HF)))))  ;
  % alpha= Nx*Ny/ ( sum(sum(real(diagHHa)))*ssigma+  sum(sum(Qp)) + sum(sum( (x-z).^2 ) ) ) 
 % alpha= Nx*Ny/ ( sum(sum( (x-z).^2 ) ) ) 
    
%Hcubic=fftshift(Hcubic);
% Nx=Nx-2*extras;
    %Ny=Ny-2*extras;

    if norm(z-zprev,'fro')^2/(Nx*Ny)<10^(-8)
        break;
    end

end


if norm(x-x_prev,'fro')^2/(Nx*Ny)<10^(-8)
    break;
end



% clear k1
% clear l1
% corrwin=5;
% fcorr=fftshift(real(ifft2(abs(fft2(gather((x))).^2))));
% fcentral=fcorr(Nx/2-corrwin:Nx/2+corrwin+1,Ny/2-corrwin:Ny/2+corrwin+1);
% %fcorrsorted=sort(fcentral(:));
% [k1, l1]=find(fcorr>mean(mean(fcentral)) & fcorr~=max(max(fcentral)));
% 
% 
% expNum=size(k1,1);
% 
% for j1=1:size(k1,1)
%     coord(j1,1)=mod(k1(j1)-Nx/2,Nx);
%    
%     coord(j1,2)=mod(l1(j1)-Ny/2,Ny);
%     
%                 rw{2*j1-1}(1)=coord(j1,1);
%                 rw{2*j1-1}(2)=coord(j1,2);    
%    % optcoord(coord(j1,1), coord(j1,2))=-1;
% end
                
disp 'End of iteration'
toc
   
end

 
%ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

%save(sprintf('file%dwithISNR%d',myinput1),'ISNR_Bayes');
 %   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (   norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )
    
 %       SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));

 7