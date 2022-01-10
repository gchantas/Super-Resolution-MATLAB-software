%function [PSNR, ISNR, SSIMres]=NonLocalPatchesFINDZ_YCbCr14_12_2018(pathHR,pathLR)


%Changes: PseudoInverse


clear all;
myinput1=1;
imNumber=4;
%randn('seed',0);



decFactor=2;
boundary=decFactor;
%pathHR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_train_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);
%pathLR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);
%pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);

 %pathHR = sprintf('/home/gchantas/Downloads/Set14x1_%d/image_%03d.png',decFactor,imNumber);
 %pathLR= sprintf('/home/gchantas/Downloads/Set14x%d/image_%03d.png',decFactor,imNumber);

pathHR = sprintf('/home/gchantas/Downloads/Set5x1_%d/image_%03d.png',decFactor,imNumber);
pathLR= sprintf('/home/gchantas/Downloads/Set5x%d/image_%03d.png',decFactor,imNumber);

g = im2double(  imread(pathLR)   );
g2=imresize(g, decFactor);

%name ='barbara.png';
%name ='Lena512.png';
%name ='boat.png';
%name =
%sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
extras=0;%3*4;%expand by $extras pixels the boundary of the image by replicating the bounds


f_= im2double(imread(pathHR));

[Nx, Ny,dummy]=size(f_);

f=f_(1:Nx,1:Ny,1);

[Nx, Ny]=size(f);


  %Nx=min(1024,Nx);
  %Ny=min(1024,Ny);
%f=f_(100:Nx+99,100:Nx+99);
f=f_(1:Nx,1:Ny);

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

    sigblurx=sigall;
    sigblury=sigall;

    for i=round(Nx/2)-12:round(Nx/2)+13
        for j=round(Ny/2)-12:round(Ny/2)+13
          % hh1(i,j)=sigblurx*cubic(-(i-Nx/2+1/4)*sigblurx-1/16)*sigblury*cubic(-(j-Ny/2+1/4)*sigblury-1/16);
                       
           %For decFactor=2
           hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);

        end
    end

hh1=hh1/sum(sum(hh1));
hh1=fftshift(hh1);
Hcubic=(fft2(hh1));

    

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
ssigma=10^(-7);

x=g2;

[x,alpha,Inv,HDDH]=StatSR(x,g,1,ssigma,Nx,Ny,decFactor,f);


%x(30:Nx-30,30:Ny-30)=(xprev(30:Nx-30,30:Ny-30));

'Stationary'
norm(f-x,'fro')^2/(Nx*Ny)


psInv=abs(Hcubic).^2;
psInv=psInv(1:Nx/2,1:Ny/2)+psInv(Nx/2+1:Nx,Ny/2+1:Ny)+psInv(1:Nx/2,Ny/2+1:Ny)+psInv(Nx/2+1:Nx,1:Ny/2);
psInv(psInv<0.000001)=0;

psInv(psInv~=0)=1./psInv(psInv~=0);
x1=real(ifft2(fft2(g).*(psInv)));
x1_=zeros(Nx,Ny);
x1_(1:decFactor:Nx,1:decFactor:Ny)=x1;
x2=decFactor^2*real(ifft2(conj(Hcubic).*fft2(x1_)));
norm(f-x2,'fro')^2/(Nx*Ny)


% 'Pseudo'
% norm(f-xPseudo,'fro')^2/(Nx*Ny)
% return
%%
Nx=Nx+2*extras;
Ny=Ny+2*extras;

hh1=zeros(Nx,Ny);
sigall=1/2;
 


%Bicubic separable filter

    sigblurx=sigall;
    sigblury=sigall;

    for i=round(Nx/2)-12:round(Nx/2)+13

        for j=round(Ny/2)-12:round(Ny/2)+13
          % hh1(i,j)=sigblurx*cubic(-(i-Nx/2+1/4)*sigblurx-1/16)*sigblury*cubic(-(j-Ny/2+1/4)*sigblury-1/16);
           %For decFactor=2
           hh1(i,j)=sigblurx*cubic(-(i-Nx/2-1+0.5)*sigblurx)*sigblury*cubic(-(j-Ny/2-1+0.5)*sigblury);
        end

    end

hh1=hh1/sum(sum(  hh1   ));
hh1=fftshift(hh1);
Hcubic=fft2(hh1);

Nx=Nx-2*extras;
Ny=Ny-2*extras;
    

%Qp=zeros(N);
   
%ISNR_Bayes=10*log10( norm(f(20:Nx-20,20:Ny-20)  - (g2(20:Nx-20,20:Ny-20)),'fro')^2/norm( f(20:Nx-20,20:Ny-20) - x1(20:Nx-20,20:Ny-20),'fro')^2  )


f_= im2double(imread(pathHR));


f=f_(1:Nx,1:Ny,1);

[Nx, Ny]=size(f);











expNum=400;
coord=zeros(expNum,2);

indx2=0;
indx=1;
optcoord=zeros(Nx,Ny);
optcoord(Nx/2+1,Ny/2+1)=1;

            for exper=1:expNum
                nu{exper}=7.0001;
                c2{exper}=12600*alpha;
                Z{exper}=ones(Nx,Ny)/expNum;
                A{exper}=zeros(Nx,Ny)/expNum;
                pof{exper}=1/expNum;
                Qp=zeros(Nx,Ny);
                E{exper}=zeros(Nx,Ny)/expNum;
            end;

            %Gaussian Blurring Function
            sigblurx=sigall*0.5;
            sigblury=sigall*0.5;
             for i=1:Nx
                for j=1:Ny
                    hh(i,j)=exp(-(abs(i-floor(Nx/2)-1).^3.0)*sigblurx).*exp(-(abs(j-floor(Ny/2)-1).^3.0)*sigblury);
                end
             end

            S{1}=1;
            hh=S{1}*(hh)/sum(sum(hh));
            HQ{1}=fft2(fftshift(hh));

            %x=g2(1:Nx,1:Ny);
            



indx2=0;
indx=1;
optcoord=zeros(Nx,Ny);
optcoord(Nx/2+1,Ny/2+1)=1;
%x=g2;
clear k1
clear l1
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
Fmask(Nx/2+1:2:Nx/2+1+50,Ny/2-26:2:Ny/2+26)=1;
Fmask=Fmask+circshift(Fmask,[1 1]);%Fmask(Nx/2+1:2:Nx/2+1+50,Ny/2-25:1:Ny/2+25)=1;


Fmask(Nx/2+1:Nx/2+3,Ny/2-3:Ny/2+4)=1;


Fmask(Nx/2+1,Ny/2+1:Ny)=0;

g3=zeros(Nx*2,Ny*2);
g3(1:Nx,1:Ny)=g2-sum(sum(g2))/(Nx*Ny);
%fcorr=fftshift(real(ifft2(abs(fft2((g2)).^2))));
fcorr=fftshift(real(ifft2(abs(fft2((g3)).^2))));
fcorr=fcorr(Nx/2+1:Nx*3/2,Ny/2+1:Ny*3/2);
%fcorrMax=fcorr(Nx+1,Ny+2);
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

% return

expNum=1;
indx=0;

for j1=1:1:200

    indx=j1;

    if mod(k1(indx)-Nx/2-1,Nx)~=0 || mod(l1(indx)-Ny/2-1,Ny)~=0
    
        coord(expNum,1)=mod(k1(indx)-Nx/2-1,Nx);
   
        coord(expNum,2)=mod(l1(indx)-Ny/2-1,Ny);

        optcoord(k1(indx), l1(indx))=-1;

        expNum=expNum+1;
    else
       [ mod(k1(indx)-Nx/2-1,Nx), mod(l1(indx)-Ny/2-1,Ny)]
    end
    
    % if  fcorr(k1(indx),l1(indx))<fcorrMax*0.93
    % break;
    % end

end

expNum=expNum-1




%for exper=1:expNum
%    E{exper}  =  gpuArray( real(  ifft2( fft2(  (g2-circshift(g2, coord(exper,:)  )).^2   )  .* conj(HQ{1})   )));
%end


   
%[loglikelihood,nu,c2]=trainPCSingleModel(expNum,x,c2,nu,coord,E,S);


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
            rw{2*exper-1}(1)=coord(exper,1);
            rw{2*exper-1}(2)=coord(exper,2);
            Z{exper}=ones(Nx,Ny)/expNum;
            A{exper}=zeros(Nx,Ny)/expNum;
            %pof{exper}=1/expNum;
            Qp=zeros(Nx,Ny);
            PW{exper}=ones(Nx,Ny)/expNum;
        end;

      % x=g2(1:Nx,1:Ny);

        ssigma=10^(-8);
    
        
     %numbers indicating by how many pixels calculation of differences goes beyond the
     %image boundary
     for exper=1:expNum
         maxX(exper)=max(abs( mod((coord(exper,1)+Nx/2),Nx)-Nx/2))+1;
         maxY(exper)=max(abs( mod((coord(exper,2)+Ny/2),Ny)-Ny/2))+1;
     end
     
    maxmaxX=max(maxX);
    maxmaxY=max(maxY);
     
     
for iter=1:10

     ZALL=zeros(Nx,Ny);

     for exper=1:expNum

        J=(S{1}/2+nu{exper}/2);

        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
      
         E   =  real(  ifft2( fft2( ( x-circshift(x, rw{2*exper-1})  ).^2  + Qp+circshift(Qp, rw{2*exper-1})) .* conj(HQ{1})  ));
      %  E   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* conj(HQ{exper})  ));

        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
        A  =   (nu{exper}+S{1})./(c2{exper}*E+nu{exper})   ;
        %Z{exper}  =   exp(  J*log(A{exper})-A{exper}.*(c2{exper}*E/2+nu{exper}/2)   );
        %Z{exper}  =   exp(  log(c2{exper})/2+(S{1}/2+nu{exper}/2-1)*(log(A{exper})-log(J)+psi(J))-A{exper}.*(c2{exper}*E/2+nu{exper}/2) -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2)  );


        Z  =   ((gamma((nu{exper}+S{1})/2)/gamma((nu{exper})/2)) *(c2{exper}/nu{exper})^0.5*(1+c2{exper}*E/nu{exper} ).^(-(nu{exper}+S{1})/2 ));%


        %c2prev{exper}=c2{exper};

        %if iter >1 
        %c2{exper}=Nx*Ny/sum(sum( (B{exper}/c2prev{exper}).*( x-circshift(x, rw{2*exper-1})  ).^2));
        %end

        %Z{exper}=(nu+c2{exper}*E).^(-J);
        B{exper}=A.*Z;
        ZALL=ZALL+Z;

     end

   % for exper=1:expNum
    %   Z{exper}=Z{exper}./ZALL;
    %   PW{exper}=(Z{exper}+.1)/(1+.1);
       % pof{exper}=sum(sum(Z{exper}));
    %end

    for exper=1:expNum
       % Z{exper}=Z{exper}./ZALL;
        %Z{exper}(find(Z{exper}<0.05))=0.0;
        B{exper}=B{exper}./ZALL;
        B{exper}=  c2{exper}  *   real(ifft2(  fft2(  B{exper}   )  .*   ( HQ{1})   ))  /   S{1};
      % B{exper}(B{exper}<100)=0;
        %[ZeroIndX(exper) ZeroIndY(exper)]=find(B{exper}<.1);
        
        rw{2*exper}  =   gpuArray2(  ssigma*B{exper}   );
    end
    
      
%     Ainv=single( real(ifft2(zeros(Nx,Ny))) );
% 
    Qp=zeros(Nx,Ny);
    Imask=ones(Nx,Ny);
    Imask(1:decFactor:Nx,1:decFactor:Ny)=0;
    Imask= not(Imask);
    hcubatrous=fftshift(real(ifft2( Hcubic )));
    hcubatrous=hcubatrous(extras+1:Nx+extras,extras+1:Ny+extras);
    hcubatrous=hcubatrous.*Imask;
   % Qp=sum(sum( real(ifft2( hcubatrous ))  .*real(ifft2( conj(hcubatrous) )) ))  /    ssigma;
    Qp= Qp+sum(sum(  hcubatrous  .*   hcubatrous))  /    ssigma;
    Qp=Qp.*Imask;
    
    for exper=1:expNum

        Qp=Qp+B{exper}+circshift(B{exper},rw{2*exper-1});%real(ifft2(conj(fft2(q1)).*(fft2(B{exper}))));
    end

    Qp=1./Qp;





    iw{1}=1;
    iw{2}=decFactor;
    iw{3}=Hcubic;
    iw{4}=extras;
    iw{5}=Nx;
    iw{6}=Ny;
    iw{7}=maxX;
    iw{8}=maxY;

    iw{9}=HDDH;

    g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
  
    x_prev=x;
    b = gpuArray2(g1);
% 
% 
%     for exper=1:expNum
%   %     rw{2*exper}=gpuArray2(padarray(rw{2*exper},[maxX(exper) maxY(exper)],'circular'));
%         ZeroIndX{exper} = find(gather(rw{2*exper})/ssigma>10);
%         %ZeroIndX(exper)=ZeroIndX(exper)+maxX(exper);
%         %ZeroIndY(exper)=ZeroIndY(exper)+maxY(exper);
%     end
%    return
%         %iw{10}=ZeroIndX;
        %iw{11}=ZeroIndY;

    % tic;   [ x, istop, itn, Anorm, Acond, rnorm, xnorm, D_r ]=cgLanczos( @Amat_fast,x, b(:), iw, rw, 0, 0, 1000, 10^(-10)); toc
     %tic;  x=PCconjGradients(x(:),@Amat_ultrafast,b(:),iw,rw,1./Qp(:), 1000,10^(-13));toc
     tic;  x=conjGradients(x(:),@Amat_ultrafast,b(:),iw,rw, 1000,10^(-8));toc

    x  =   reshape(x,  Nx, Ny );

    % D_r=reshape(D_r,  Nx, Ny );


    %ISNR_Bayes=10*log10(norm(fcolor(20:Nx-20,20:Ny-20)-g__2(20:Nx-20,20:Ny-20),'fro')^2/(norm(fcolor(20:Nx-20,20:Ny-20,1)-x_r(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,2)-x_g(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,3)-x_b(20:Nx-20,20:Ny-20),'fro')^2))

    boundary=decFactor+1;
    
    %ssigma=(norm(g-myimresize(x,Nx,Ny,Hcubic,decFactor,extras),'fro')^2)/(Nx*Ny/decFactor^2)
    ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))
    PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - double(uint8(255*x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)))/255,   'fro'    )^2 ) )

    SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));

iter
if norm(x-x_prev,'fro')<0.0001
    return;
end

end
 
%ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

%save(sprintf('file%dwithISNR%d',myinput1),'ISNR_Bayes');
 %   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (   norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )
    
 %       SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));

