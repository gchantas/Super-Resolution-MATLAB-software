function [PSNR, ISNR, SSIMres]=SRviaDenoisingYCbCrGauss24_1_2020(pathHR,pathLR,decFactor)
%clear all;
%myinput1=1;
%imNumber=13;
%randn('seed',0);



%decFactor=2;

%pathHR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_train_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);
%pathLR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);

%pathHR = sprintf('/home/gchantas/Downloads/Set14x1_%d/image_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Downloads/Set14x%d/image_%03d.png',decFactor,imNumber);



%name ='barbara.png';
%name ='Lena512.png';
%name ='boat.png';
%name = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
extras=0;%3*4;%expand by $extras pixels the boundary of the image by replicating the bounds


f_= im2double(imread(pathHR));

[Nx, Ny,dummy]=size(f_);

f=f_(1:Nx,1:Ny,1);

[Nx, Ny]=size(f);


  %Nx=min(1024,Nx);
  %Ny=min(1024,Ny);
%f=f_(100:Nx+99,100:Nx+99);
f=f_(1:Nx,1:Ny);

g = im2double(imread(pathLR));


Nx=Nx+2*extras;
Ny=Ny+2*extras;

hh1=zeros(Nx,Ny);
sigall=1/decFactor;
 

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

Nx=Nx-2*extras;
Ny=Ny-2*extras;
    

norm(myimresize(f,Nx,Ny,Hcubic,decFactor,extras)-imresize(f,1/decFactor),'fro')^2/(Nx*Ny/(decFactor*2))




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

g2=imresize(g,decFactor);

alpha=.1;
ssigma=.00001;

dcf=decFactor;
decFactor=decFactor^2/2;

I=makeBlockDiagMatCells(Nx+2*extras,Ny+2*extras,decFactor*2,0);%%Block diagonal identity matrix
HcubBlock=makeDiagMatCells(Hcubic,Nx+2*extras,Ny+2*extras,decFactor*2,0);
HcubBlockConj=makeDiagMatCells(conj(Hcubic),Nx+2*extras,Ny+2*extras,decFactor*2,0);

C=multiplyBlockDiagMat(HcubBlockConj,multiplyBlockDiagMat(I,HcubBlock,decFactor*2),decFactor*2);


g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,dcf,extras) ;

G_=decFactor*padarray(g1,[extras extras],'symmetric');

G_=fft2(G_);

z_=decFactor*padarray(g1,[extras extras],'symmetric');

Z_=fft2(z_);

%g3=imresize(g,2);


%tic

Ia=makeDiagMatCells(  alpha*ssigma*abs(Q).^2,Nx+2*extras,Ny+2*extras,decFactor*2,0   );%%diagonal identity matrix (no blocks)
C2=addBlockMat(C,Ia,decFactor*2);
Inv=invertBlockMatrix(C2,decFactor*2);

G3=makeBlockVec( G_ , Nx+2*extras,Ny+2*extras,decFactor,0);%applyBlockDiag(HcubBlockConj,G2_,decFactor*2);
bVec=makeBlockVec(  Z_ , Nx+2*extras,Ny+2*extras,decFactor,0);


xblock=applyBlockDiag(Inv,bVec,decFactor*2);

vec=deblockVec(xblock,Nx+2*extras,Ny+2*extras,decFactor*2,0);

x=real(ifft2(vec));
x=x(extras+1:Nx+extras,extras+1:Ny+extras);
%diagHHa=tageDiagBlockMat(Inv,Nx,Ny,decFactor*2,0);
xprev=x;
decFactor=dcf;
%boundaryX=20;
%x(boundaryX:Nx-boundaryX,boundaryX:Ny-boundaryX)=gather(xprev(boundaryX:Nx-boundaryX,boundaryX:Ny-boundaryX));

x=g2;

% if decFactor==2
%     ssigma=10^(-7);
% [x,alpha1,ssigma,Inv,HDDH]=StatSR(x,g,1,ssigma,Nx,Ny,decFactor,f);
% elseif decFactor==4
%     ssigma=10^(-6);
%     [x,alpha1,ssigma,Inv,HDDH]=StatSRx4(x,g,1,ssigma,Nx,Ny,decFactor,f);
% %alpha=alpha*2;
% end


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
%name = sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);


%g=g_(60/2:min([Nx Ny])/2+60/2-1,60/2:min([Nx Ny])/2+60/2-1);
%g=g_;
%%%%%%%%%%%%%%%%%%%%%%
%PSNRest=10*log10( 1/(ssigma*255) )

  
    
%[fhat,a,ssigma]=EM(single(g),single(real(ifft2(H))),single(real(ifft2(Q))),N,N,100);
%[a,ssigma,fhat]=stat_rest(1,1,g,Q,H,N);
% a=50000;
%ssigma



%Mfg=conj(H).*fft2(g)./(Hf+a*ssigma*Qf/10);
%fhat=real(ifft2(Mfg));

%isnrstat=10*log10(norm(f-g,'fro')^2/norm(f-fhat,'fro')^2)

%PSNRtrue=10*log10( NN/norm(g-f,'fro')^2 )

%Qp=zeros(N);
%ISNR_Bayes=10*log10( norm(f(20:Nx-20,20:Ny-20)  - (g2(20:Nx-20,20:Ny-20)),'fro')^2/norm( f(20:Nx-20,20:Ny-20) - x1(20:Nx-20,20:Ny-20),'fro')^2  )



f_= im2double(imread(pathHR));

f=f_(1:Nx,1:Ny,1);
[Nx, Ny]=size(f);

g = im2double(imread(pathLR));


expNum=100;
coord=zeros(expNum,2);



indx2=0;
indx=1;
optcoord=zeros(Nx,Ny);
optcoord(Nx/2+1,Ny/2+1)=1;
%x=g2;
clear k1
clear l1
Fmask=zeros(Nx,Ny);
Fmask(Nx/2+1:2:Nx/2+1+90,Ny/2-46:2:Ny/2+46)=1;
Fmask=Fmask+circshift(Fmask,[1 1]);%Fmask(Nx/2+1:2:Nx/2+1+50,Ny/2-25:1:Ny/2+25)=1;


Fmask(Nx/2+1:Nx/2+8,Ny/2-4:Ny/2+4)=1;


Fmask(Nx/2+1,Ny/2+1:Ny)=0;

g3=zeros(Nx*2,Ny*2);
g3(1:Nx,1:Ny)=g2;
%g3(1:Nx,1:Ny)=g3(1:Nx,1:Ny);
%fcorr=fftshift(real(ifft2(abs(fft2((g2)).^2))));
fcorr=fftshift(real(ifft2(abs(fft2((g3-mean(mean(g3)))).^2))));
fcorr=fcorr(Nx/2+1:Nx*3/2,Ny/2+1:Ny*3/2);
fcorrMax=fcorr(Nx/2+1,Ny/2+1);

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

    if mod(k1(indx)-Nx/2-1,Nx)~=0 || mod(l1(indx)-Ny/2-1,Ny)~=0
    
        coord(expNum,1)=mod(k1(indx)-Nx/2-1,Nx);
   
        coord(expNum,2)=mod(l1(indx)-Ny/2-1,Ny);

        optcoord(k1(indx), l1(indx))=-1;

        expNum=expNum+1;
    else
       [ mod(k1(indx)-Nx/2-1,Nx), mod(l1(indx)-Ny/2-1,Ny)]
    end
    
    if  fcorr(k1(indx),l1(indx))<fcorrMax*0.60
        break;
    end

end

expNum=expNum-1


%indx2=0;
%indx=1;
%for j1=size(fcorrsorted,1)-(corrwin+1)^2:size(fcorrsorted,1)
%    ind1=fcorr(fcorr==fcorrsorted(j1));
%    [k1, l1]=ind2sub(size(f), ind1)
    
    
%end
%  
% figure,imagesc(fcorr);
% figure,imagesc(fcorrsorted);

% for j1=1:1:5
%     for k1=-j1:round(abs(j1)^.4):j1+round(abs(j1)^.4)-1
%         for l1=-j1:round(abs(j1)^.4):j1+round(abs(j1)^.4)-1
%         
%   %            for k1=-j1:j1
%   %              for l1=-j1:j1
%              % if(k1~=0&&l1~=0)
%             %(k1>0 | (k1==0 & l1>0))&
%             
%             if (  k1>0  ||   (k1==0 && l1<0))  &&   (abs(k1)  ==   j1  ||   abs(l1)  ==   j1   )%(k1~=0 | l1~=0)%sqrt(k1^2+l1^2)<6&
%             %   if (abs(k1)==j1 | abs(l1)==j1)&(k1~=0 | l1~=0)
%                 % l=l1
%                 %dirDiff(indx)  =   (norm( (W.*abs( fcolor(:,:,1)-circshift(fcolor(:,:,1), [k1 l1]  )   )).^2, 'fro' )  +   norm( (W.*abs( fcolor(:,:,2)-circshift(fcolor(:,:,2), [k1 l1]) )  ).^2, 'fro' )  +   norm( (W.*abs( fcolor(:,:,3)-circshift(fcolor(:,:,3), [k1 l1])  ) ).^2, 'fro' ));
% 
%             
%             
%                 sigma1 = 1;
%                 sigma2 = 1;
%                 scale1 = 1;
%                 scale2 = 1;
%                 sigma1 = scale1*sigma1;
%                 sigma2 = scale2*sigma2;
%                 coord(indx,1)=mod(k1,Nx);
%                 coord(indx,2)=mod(l1,Ny);
%                 
%                 Theta = atan(l1/k1);
%  
%                 hh=zeros(Nx,Ny);
% 
%                 sigall=1;
% 
%                
%                 %hh=imrotate(hh,Theta*180/(pi),'crop');
%                 %figure,imagesc(hh);
% 
%               
%   
%                 %if k1==0 && l1==0
%                 %optcoord(l1+Nx/2+1,k1+Ny/2+1)=-1;
% 
%                 %end
%                     
%                     
%             indx=indx+1;
%             end
% 
%         end
%     end
% end
sigall=1;
%sigall=sqrt(abs(k1^2+l1^2));
            %Gaussian Blurring Function
            sigblurx=sigall*0.25;
            sigblury=sigall*0.25;
             for i=1:Nx
                for j=1:Ny
                    hh(i,j)=exp(  -(abs(i-floor(Nx/2)-1).^3.0   )*sigblurx).*exp(  -(abs(j-floor(Ny/2)-1).^3.0)*sigblury   );
                end
             end


            S{1}=1;
            hh=S{1}*(hh)/sum(sum(hh));
            HQ{1}=gpuArray2(fft2(fftshift(hh)));


           for exper=1:expNum
                nu{exper}=4.0001;
                c2{exper}=5000.00;
                rw{2*exper-1}(1)=coord(exper,1);
                rw{2*exper-1}(2)=coord(exper,2);
                Z{exper}=gpuArray2(ones(Nx,Ny)/expNum);

                A=gpuArray2(zeros(Nx,Ny)/expNum);
                pof{exper}=1/expNum;
                Qp=gpuArray2(zeros(Nx,Ny));
                                  
             %   E{exper} =  gpuArray2(real(  ifft2( fft2( (f-circshift(f, rw{2*exper-1})  ).^2 ) .* conj(HQ{1}))  ));

                %PW{exper}=ones(Nx,Ny)/expNum;
            end;
            
            
            %  [loglikelihood,nu,c2]=trainPCSingleModel(expNum,x,c2,nu,coord,E,S);% 
            %  [loglikelihood,nu,c2]=trainPCMixture(expNum,x,c2,nu,coord,E,S);% 

          %  expNum=indx-1;
%return


E=zeros(Nx,Ny);

%g2_=imsharpen(g2,'Radius',3,'Amount' ,3);

z=gpuArray2(x);

%x=g2;



alpha=10;
ssigma=.0000001;

dcf=decFactor;
decFactor=decFactor^2/2;

I=makeBlockDiagMatCells(Nx+2*extras,Ny+2*extras,decFactor*2,0);%%Block diagonal identity matrix

HcubBlock=makeDiagMatCells(Hcubic,Nx+2*extras,Ny+2*extras,decFactor*2,0);
HcubBlockConj=makeDiagMatCells(conj(Hcubic),Nx+2*extras,Ny+2*extras,decFactor*2,0);

%G2_=subfoldVec(G_,decFactor*2,0);
C=multiplyBlockDiagMat(HcubBlockConj,multiplyBlockDiagMat(I,HcubBlock,decFactor*2),decFactor*2);

g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,dcf,extras);
G_=decFactor*padarray(g1,[extras extras],'symmetric')*2;

G_=fft2(G_);
Imask=ones(Nx,Ny);
Imask(1:decFactor:Nx,1:decFactor:Ny)=0;
Imask= not(Imask);
decFactor=dcf;


for iter=1:200
tic
    
     ZALL=zeros(Nx,Ny);
%tic
     for exper=1:expNum

        J=(S{1}/2+nu{exper}/2);

%         qa=zeros(Nx,Ny);
%         qa(1,1)=1;
%         qa(coord(exper,1)+1,coord(exper,2)+1)=1;
% 
%         q1=real(ifft2(fft2(qa))).^2;
        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
      % Qe =  gpuArray2(real(  ifft2( fft2( ( z-circshift(z, rw{2*exper-1})  ))))).^2;

          E =  gpuArray2(real(  ifft2( fft2( ( x-circshift(x, rw{2*exper-1})  ).^2  + Qp+circshift(Qp, rw{2*exper-1})) .* conj(HQ{1})  )));
  
      %  E   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* conj(HQ{exper})  ));
       
        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
           
%         c2prev{exper}=c2{exper};
%         
%         if iter >1
%              c2{exper}=Nx*Ny/sum(sum( (B{exper}).*E));
%         end
       

     % A=real(ifft2(  fft2(  (nu+S{1})./(nu+c2{exper}*Qe) )  .*   conj(HQ{1}) ));
     % A =   (nu+S{1})./(c2{exper}*E+nu);

     % Z{exper} = real(ifft2(  fft2(  (nu+c2{exper}*Qe).^(-J) )  .*   conj(HQ{1}) ));
     % Z{exper}  =    exp(  J*log(A)-A.*(c2{exper}*E/2+nu/2)   );
     % Z{exper}  = (c2{exper})^0.5*(1+c2{exper}*E/nu ).^(-(nu+S{1})/2 );%

    A =  ones(Nx,Ny);%(nu{exper}+S{1})./(c2{exper}*E+nu{exper});% 
        %Z{exper}  =   exp(  J*log(A{exper})-A{exper}.*(c2{exper}*E/2+nu{exper}/2)   );
      %  Z{exper}  =   exp(  log(c2{exper})/2+(S{1}/2+nu{exper}/2-1)*(log(A{exper})-log(J)+psi(J))-A{exper}.*(c2{exper}*E/2+nu{exper}/2) -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2)  );
     % Z{exper}  =  ( (gamma((nu{exper}+S{1})/2)/gamma((nu{exper})/2)) *(c2{exper}/nu{exper})^0.5*(1+c2{exper}*E/nu{exper} ).^(-(nu{exper}+S{1})/2 ));%
      Z{exper}=exp(-0.5*c2{exper}*E);


      %  Z{exper}=(nu+c2{exper}*E).^(-J);

        %B{exper}=A{exper};
        ZALL=ZALL+Z{exper};
       %B{exper}=A.*Z{exper};
       B{exper}=A;
    end

    for exper=1:expNum
        Z{exper}=Z{exper}./ZALL;
   %  PW{exper}=(Z{exper}+.002)/(1+.002);

        %Z{exper}(find(Z{exper}<0.05))=0.0;
        B{exper}=B{exper}.*Z{exper};
       %B{exper}=B{exper}.*(Z{exper}+0.01)/(1+0.01);
       % c2{exper}=Nx*Ny/sum(sum(B{exper}.*E{exper}));

        %c2{exper}=Nx*Ny/( sum(sum(B{exper}.*( x-circshift(x, rw{2*exper-1})  ).^2))+ );
        B{exper}=   gpuArray2( c2{exper}* real(ifft2(  fft2(  B{exper}   )  .*   ( HQ{1})   ))  /   S{1});
        
        %B{exper}=B{exper}.*not(Imask);
       % rw{2*exper}  =   gpuArray(  B{exper}   );
    end
        %disp 'time to B, E and Z'
%toc
    
%x=g2;

alpha=2000*sqrt(iter);%.^1.1;

ssigma=10^(-7);

x_prev=x;

for iter2=1:100

    iw{1}=alpha;
    iw{2}=decFactor;
    iw{3}=ones(Nx,Ny);
    iw{4}=extras;
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


   
       %PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - z(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )
    %          if norm(z-z1,'fro')<0.0001
    %              'Number of iteration vardenoise stopped'
    %              iter3
    %              break
    %          end
    
       % z=gpuArray2(g2);

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
        B{exper}=padarray(B{exper},[maxX maxY],'circular');
    end

    %tic

for iter3=1:1
zprev=z;

        [z,Ball]=VariationalDenoise(padarray(z1,[maxX maxY],'circular'),padarray(x,[maxX maxY],'circular'),B,Nx,Ny,coord,expNum,alpha,maxX,maxY);

        if norm(z-z1,'fro')^2/(Nx*Ny)<10^(-7)
            norm(z-z1,'fro')
            disp('z evaluation converged')
            break;
        end;
        z1=z;
         
end



    for exper=1:expNum
        B{exper}=B{exper}(maxX+1:Nx+maxX,maxY+1:Ny+maxY);
    end
 
   % toc

    'time to var x'
   
    Qp=1./(alpha+Ball);
    
   
   %PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - z(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )

   
   %z=reshape(z,  Nx, Ny );
%    g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
 
 %  b = gpuArray2(g1+alpha*ssigma*z);
   
  % x1=real(ifft2( HF.* padarray(b,[extras extras],'circular')));
  % x=x1(extras+1:Nx-extras,extras+1:Ny-extras);

  iw{1}=alpha*ssigma;
  iw{2}=decFactor;
  iw{3}=Hcubic;
  iw{4}=extras;

%     g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
%     b = gpuArray2(g1+alpha*ssigma*z);
%  tic
%     xprev=x;
%     x=cgLanczos( @EasySRAmat,gpuArray2(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 1000, 0.000000001); 
%  toc
%     x=reshape(x,  Nx, Ny );
%alpha=1;
z_=padarray(z,[extras extras],'symmetric');
Z_=fft2(z_);
%g3=imresize(g,2);

%tic

decFactor=decFactor^2/2;

Ia=makeDiagMatCells(  alpha*ssigma*ones(Nx+2*extras,Ny+2*extras),Nx+2*extras,Ny+2*extras,decFactor*2,0   );%%diagonal identity matrix (no blocks)
C2=addBlockMat(C,Ia,decFactor*2);
Inv=invertBlockMatrix(C2,decFactor*2);

G3=makeBlockVec( G_ ,Nx+2*extras,Ny+2*extras,decFactor,0);%applyBlockDiag(HcubBlockConj,G2_,decFactor*2);
bVec=makeBlockVec(  ssigma*alpha*Z_,Nx+2*extras,Ny+2*extras,decFactor,0);
bVec2=addBlockVec(G3,bVec,decFactor*2);

xblock=applyBlockDiag(Inv,bVec2,decFactor*2);

vec=deblockVec(xblock,Nx+2*extras,Ny+2*extras,decFactor*2,0);
xprev=x;
x=real(ifft2(vec));
x=x(extras+1:Nx+extras,extras+1:Ny+extras);

decFactor=dcf;
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
%     
%     
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

