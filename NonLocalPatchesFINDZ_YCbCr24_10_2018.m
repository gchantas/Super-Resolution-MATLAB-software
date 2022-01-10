%function [PSNR, ISNR, SSIMres]=NonLocalPatchesYCbCr24_10_2018(imNumber)
clear all;
myinput1=1;
imNumber=1;
%randn('seed',0);

experimentN=3;


decFactor=2;
boundary=decFactor;
%pathHR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_train_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);
%pathLR = sprintf('/home/gchantas/Documents/SR/DIV2Kdataset/DIV2K_valid_LR_bicubic/X2/%04dx%d.png',imNumber,decFactor);
%pathHR = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
%pathLR= sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);

pathHR = sprintf('/home/gchantas/Downloads/Set14x1_%d/image_%03d.png',decFactor,imNumber);
pathLR= sprintf('/home/gchantas/Downloads/Set14x%d/image_%03d.png',decFactor,imNumber);


%name ='barbara.png';
%name ='Lena512.png';
%name ='boat.png';
%name = sprintf('/home/gchantas/Documents/SR/Urban100x1_%d/img_%03d.png',decFactor,imNumber);
extras=3*4;%expand by $extras pixels the boundary of the image by replicating the bounds


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
        
%%%Generate BLURING to produce observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

hh1=hh1/sum(sum(hh1));
hh1=fftshift(hh1);
Hcubic=conj(fft2(hh1));

Nx=Nx-2*extras;
Ny=Ny-2*extras;
    

norm(myimresize(f,Nx,Ny,Hcubic,decFactor,extras)-imresize(f,1/decFactor),'fro')^2/(Nx*Ny)


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
%name = sprintf('/home/gchantas/Documents/SR/Urban100x%d/img_%03d.png',decFactor,imNumber);

g = im2double(imread(pathLR));

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




f_= im2double(imread(pathHR));


f=f_(1:Nx,1:Ny,1);

[Nx, Ny]=size(f);



g = im2double(imread(pathLR));



expNum=400;
coord=zeros(expNum,2);

indx2=0;
indx=1;
optcoord=zeros(Nx,Ny);
optcoord(Nx/2+1,Ny/2+1)=1;

            for exper=1:expNum
                nu{exper}=2.0001;
                c2{exper}=1000.01;

                Z{exper}=ones(Nx,Ny)/expNum;

                A{exper}=zeros(Nx,Ny)/expNum;
                pof{exper}=1/expNum;
                Qp=zeros(Nx,Ny);
                E{exper}=zeros(Nx,Ny)/expNum;

            end;
            
            
            %Gaussian Blurring Function
            sigblurx=sigall*0.4;
            sigblury=sigall*0.4;
             for i=1:Nx
                for j=1:Ny
                    hh(i,j)=exp(-(abs(i-floor(Nx/2)-1).^3.0)*sigblurx).*exp(-(abs(j-floor(Ny/2)-1).^3.0)*sigblury);
                end
             end


            S{1}=1;
            hh=S{1}*(hh)/sum(sum(hh));
            HQ{1}=fft2(fftshift(hh));
            

    
            x=g2(1:Nx,1:Ny);
         
            loglikelihood=-inf;
            
            ZALL=zeros(Nx,Ny);

loglikelihood_vec(1)=-inf;





for j1=1:2:3
    for k1=-j1:round(abs(j1)^.4):j1+round(abs(j1)^.4)-1
        for l1=-j1:round(abs(j1)^.4):j1+round(abs(j1)^.4)-1
        
  %            for k1=-j1:j1
  %              for l1=-j1:j1
             % if(k1~=0&&l1~=0)
            %(k1>0 | (k1==0 & l1>0))&
            
            if (  k1>0  ||   (k1==0 && l1<0))  &&   (abs(k1)  ==   j1  ||   abs(l1)  ==   j1   )%(k1~=0 | l1~=0)%sqrt(k1^2+l1^2)<6&
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

        %hh=imrotate(hh,Theta*180/(pi),'crop');
        %figure,imagesc(hh);
        %if k1==0 && l1==0
        %end

        indx=indx+1;
        expNum=indx-1;

        E{expNum}  =  gpuArray( real(  ifft2( fft2(  (f-circshift(f, coord(expNum,:)  )).^2   )  .* conj(HQ{1})   )));

        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
        %E   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* conj(HQ{exper})  ));

        %for emiter=1:10
      %  prevlikelihood=loglikelihood;
%        [loglikelihood,nu,c2]=trainPCMixture(expNum,x,c2,nu,coord,E,S);
%        % disp('Likelihoods')
%        % prevlikelihood
%        % loglikelihood
%  disp('Nus')
% 
%          nu{1:expNum}

           %if prevlikelihood>loglikelihood
%            if nu{expNum}<.1 %|| prevlikelihood>loglikelihood
%              %loglikelihood  =   prevlikelihood ;
%              indx=indx-1;
%              disp('nu excluded')
%              nu{expNum}
%   
%            else
%             % loglikelihood_vec(expNum)=gather(loglikelihood);
              optcoord(k1+Nx/2+1,l1+Ny/2+1)=-1;
% 
%          end

        end
    end
    end
end


        expNum=indx-1;
        [loglikelihood,nu,c2]=trainPCMixture(expNum,x,c2,nu,coord,E,S);

        disp('Nus')

        nu{1:expNum}

return;
        %g2_=imsharpen(g2,'Radius',3,'Amount' ,3);

        for exper=1:expNum
            %nu=0.00001;
            %c2{exper}=1.01;
            rw{2*exper-1}(1)=coord(exper,1);
            rw{2*exper-1}(2)=coord(exper,2);
            Z{exper}=ones(Nx,Ny)/expNum;

            A{exper}=zeros(Nx,Ny)/expNum;
            pof{exper}=1/expNum;
            Qp=zeros(Nx,Ny);

        end;


        x=g2(1:Nx,1:Ny);



ssigma=10^(-7);

for iter=1:10

     ZALL=zeros(Nx,Ny);

     for exper=1:expNum

        J=(S{1}/2+nu{exper}/2);

        qa=zeros(Nx,Ny);
        qa(1,1)=1;
        qa(coord(exper,1)+1,coord(exper,2)+1)=1;

        q1=real(ifft2(fft2(qa))).^2;
        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
      
         E   =  real(  ifft2( fft2( ( x-circshift(x, rw{2*exper-1})  ).^2  + Qp+circshift(Qp, rw{2*exper-1})) .* conj(HQ{1})  ));
      %  E   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* conj(HQ{exper})  ));

        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
        A{exper}  =   (nu{exper}+S{1})./(c2{exper}*E+nu{exper})   ;
        %Z{exper}  =   exp(  J*log(A{exper})-A{exper}.*(c2{exper}*E/2+nu{exper}/2)   );
        Z{exper}  =   exp(  log(c2{exper})/2+(S{1}/2+nu{exper}/2-1)*(log(A{exper})-log(J)+psi(J))-A{exper}.*(c2{exper}*E/2+nu{exper}/2) -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2)  );


% 
%         c2prev{exper}=c2{exper};
% 
%         if iter >1 
%              c2{exper}=Nx*Ny/sum(sum( (B{exper}/c2prev{exper}).*( x-circshift(x, rw{2*exper-1})  ).^2));
%         end

        %Z{exper}=(nu+c2{exper}*E).^(-J);
        %B{exper}=A{exper};
        ZALL=ZALL+Z{exper};

    end

    for exper=1:expNum
        Z{exper}=Z{exper}./ZALL;
        %Z{exper}(find(Z{exper}<0.05))=0.0;
        B{exper}=A{exper}.*Z{exper};
        B{exper}=  c2{exper}  *   real(ifft2(  fft2(  B{exper}   )  .*   ( HQ{1})   ))  /   S{1};
        rw{2*exper}  =   gpuArray(  ssigma*B{exper}   );
    end
    


    Ainv=single( real(ifft2(zeros(Nx,Ny))) );

    hcubatrous=real(ifft2( Hcubic ));
    hcubatrous(1:decFactor:Nx,1:decFactor:Ny)=0;
    %Qp=sum(sum( real(ifft2( hcubatrous ))  .*real(ifft2( conj(hcubatrous) )) ))  /    ssigma;
    Qp= sum(sum( hcubatrous .*hcubatrous))  /    ssigma;

    for exper=1:expNum


        qa=zeros(Nx,Ny);
        qa(1,1)=1;

        qa(coord(exper,1)+1,coord(exper,2)+1)=1;
        q1=real(ifft2((fft2(qa)))).^2;
        Qp=Qp+real(ifft2(conj(fft2(q1)).*(fft2(B{exper}))));
    end

    Qp=1./Qp;


    itnum=0;


    maxX=40;
    maxY=40;

    iw{1}=1;
    iw{2}=decFactor;
    iw{3}=Hcubic;
    iw{4}=extras;
    iw{5}=Nx;
    iw{6}=Ny;
    iw{7}=maxX;
    iw{8}=maxY;

    g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
  
    x_prev=x;
    b = gpuArray(g1);
    
      
    for exper=1:expNum
        rw{2*exper}=gpuArray(padarray(rw{2*exper},[maxX maxY],'circular'));
    end
    
    
    tic;   [ x, istop, itn, Anorm, Acond, rnorm, xnorm, D_r ]=cgLanczos( @Amat_fast,gpuArray(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 1200, 10^(-16)); toc

    x=reshape(x,  Nx, Ny );
   % D_r=reshape(D_r,  Nx, Ny );


    %ISNR_Bayes=10*log10(norm(fcolor(20:Nx-20,20:Ny-20)-g__2(20:Nx-20,20:Ny-20),'fro')^2/(norm(fcolor(20:Nx-20,20:Ny-20,1)-x_r(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,2)-x_g(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,3)-x_b(20:Nx-20,20:Ny-20),'fro')^2))

    boundary=decFactor;
    

%ssigma=(norm(g-myimresize(x,Nx,Ny,Hcubic,decFactor,extras),'fro')^2)/(Nx*Ny/decFactor^2)


    ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

    PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )


    SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));

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
iter
if norm(x-x_prev,'fro')<0.0001
    return;
end

end
 
%ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

%save(sprintf('file%dwithISNR%d',myinput1),'ISNR_Bayes');
 %   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (   norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )
    
 %       SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));

