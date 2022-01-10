%function [PSNR, ISNR, SSIMres]=SRviaDenoisingYCbCr31_10_2018(imNumber)
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



g2=imresize(g,decFactor);


%ISNR_Bayes=10*log10( norm(f(20:Nx-20,20:Ny-20)  - (g2(20:Nx-20,20:Ny-20)),'fro')^2/norm( f(20:Nx-20,20:Ny-20) - x1(20:Nx-20,20:Ny-20),'fro')^2  )




f_= im2double(imread(pathHR));


f=f_(1:Nx,1:Ny,1);

[Nx, Ny]=size(f);



g = im2double(imread(pathLR));



expNum=100;
coord=zeros(expNum,2);

indx2=0;
indx=1;
%optcoord=zeros(Nx,Ny);
%optcoord(Nx/2+1,Ny/2+1)=1;

for j1=1:2:5
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
                %optcoord(l1+Nx/2+1,k1+Ny/2+1)=-1;

                %end
                    
                    
            indx=indx+1;
            end

        end
    end
end

            %Gaussian Blurring Function
            sigblurx=sigall*0.4;
            sigblury=sigall*0.4;
             for i=1:Nx
                for j=1:Ny
                    hh(i,j)=exp(  -(abs(i-floor(Nx/2)-1).^3.0   )*sigblurx).*exp(  -(abs(j-floor(Ny/2)-1).^3.0)*sigblury   );
                end
             end


            S{1}=7;
            hh=S{1}*(hh)/sum(sum(hh));
            HQ{1}=gpuArray(fft2(fftshift(hh)));

            expNum=indx-1;

            for exper=1:expNum
                nu=0.0001;
                c2{exper}=1.000;
                rw{2*exper-1}(1)=coord(exper,1);
                rw{2*exper-1}(2)=coord(exper,2);
                Z{exper}=gpuArray(ones(Nx,Ny)/expNum);

                A=gpuArray(zeros(Nx,Ny)/expNum);
                pof{exper}=1/expNum;
                Qp=gpuArray(zeros(Nx,Ny));
                          

            end;
E=zeros(Nx,Ny);

%g2_=imsharpen(g2,'Radius',3,'Amount' ,3);
x=gpuArray(g2);

z=gpuArray(g2);

alpha=100;
ssigma=0.0000001;
tic
for iter=1:250
tic
    
     ZALL=zeros(Nx,Ny);

     for exper=1:expNum

        J=(S{1}/2+nu/2);

%         qa=zeros(Nx,Ny);
%         qa(1,1)=1;
%         qa(coord(exper,1)+1,coord(exper,2)+1)=1;
% 
%         q1=real(ifft2(fft2(qa))).^2;
        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
      
         E =  gpuArray(real(  ifft2( fft2( ( z-circshift(z, rw{2*exper-1})  ).^2  + Qp+circshift(Qp, rw{2*exper-1})) .* conj(HQ{1})  )));
      %  E   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* conj(HQ{exper})  ));
       
        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
           
%         c2prev{exper}=c2{exper};
%         
%         if iter >1 
%              c2{exper}=Nx*Ny/sum(sum( (B{exper}).*E));
%         end     
        A  =   (nu+S{1})./(c2{exper}*E+nu)   ;
        %Z{exper}  =   exp(  J*log(A)-A.*(c2{exper}*E/2+nu/2)   );


        
        
        Z{exper}=(nu+c2{exper}*E).^(-J);

        %B{exper}=A{exper};
        ZALL=ZALL+Z{exper};
       B{exper}=A.*Z{exper};
    end

    for exper=1:expNum
        %Z{exper}=Z{exper}./ZALL;
        %Z{exper}(find(Z{exper}<0.05))=0.0;
        B{exper}=B{exper}./ZALL;
       % c2{exper}=Nx*Ny/sum(sum(B{exper}.*E{exper}));

        %c2{exper}=Nx*Ny/( sum(sum(B{exper}.*( x-circshift(x, rw{2*exper-1})  ).^2))+ );
        B{exper}=   gpuArray( c2{exper}* real(ifft2(  fft2(  B{exper}   )  .*   ( HQ{1})   ))  /   S{1});
        rw{2*exper}  =   gpuArray(  B{exper}   );
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
% 
%     Ainv=single( real(ifft2(zeros(Nx,Ny))) );
% %    xprev=x;
%     
%     Qp=sum(sum( real(ifft2( Hcubic ))  .*real(ifft2( conj(Hcubic) )) ))  /    ssigma;
% 
%     for exper=1:expNum
% 
%         %pof{exper}=exp( psi(Z{exper}+10000) -psi(A0) );%(Z{exper}+0.0812)./ A0;
%         qa=zeros(Nx,Ny);
%         qa(1,1)=1;
% 
%         qa(coord(exper,1)+1,coord(exper,2)+1)=1;
%         q1=real(ifft2((fft2(qa)))).^2;
%         Qp=Qp+real(ifft2(conj(fft2(q1)).*(fft2(B{exper}))));
%     end
% 
%     Qp=1./Qp;

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



alpha=100*sqrt(iter);%.^1.1;
ssigma=0.0000001;

for iter2=1:1
    
    iw{1}=alpha;
    iw{2}=decFactor;
    iw{3}=ones(Nx,Ny);
    iw{4}=extras;
    %imresize(  g,   decFactor);
    %g1(2:decFactor:Nx, 2:decFactor:Ny)=g;

    %g1=zeros(Nx,Ny);
    %g1(1:decFactor:Nx,1:decFactor:Ny)=g_r;
    %g1 = real(ifft2( fft2(g1) .* Hcubic));%imresize(g_r,decFactor);
   % g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
  
    x_prev=x;
    b = gpuArray(x)*alpha;
    %tic;   [ z, istop, itn, Anorm, Acond, rnorm, xnorm, D_r ]=cgLanczos( @DenoiseAmat,gpuArray(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 100, 0.0000001); toc




    z=gpuArray(g2);
    z1=gpuArray(z);

    maxX=20;
    maxY=20;
    for exper=1:expNum
        B{exper}=padarray(B{exper},[maxX maxY],'circular');
    end
tic
    for iter3=1:10

        
        [z,Ball]=VariationalDenoise(padarray(z1,[maxX maxY],'circular'),padarray(x,[maxX maxY],'circular'),B,Nx,Ny,coord,expNum,alpha,maxX,maxY);

    %PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - z(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )
    %          if norm(z-z1,'fro')<0.0001
    %              'Number of iteration vardenoise stopped'
    %              iter3
    %              break
    %          end

    z1=z;
         
    end

    for exper=1:expNum
        B{exper}=B{exper}(maxX+1:Nx+maxX,maxY+1:Ny+maxY);
    end
 
    toc


    'time to var x'
    
 

   
   Qp=1./(alpha+Ball);
    
   
   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - z(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )

   
   %z=reshape(z,  Nx, Ny );
%    g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
 
 %  b = gpuArray(g1+alpha*ssigma*z);
   
  % x1=real(ifft2( HF.* padarray(b,[extras extras],'circular')));
  % x=x1(extras+1:Nx-extras,extras+1:Ny-extras);
   iw{1}=alpha*ssigma;
   iw{2}=decFactor;
    iw{3}=Hcubic;
    iw{4}=extras;
    g1 = myimresizeTranspose(g,Nx,Ny,Hcubic,decFactor,extras) ;
 
    b = gpuArray(g1+alpha*ssigma*z);
 
    x=cgLanczos( @EasySRAmat,gpuArray(zeros(Nx,Ny)), b(:), iw, rw, 0, 0, 100, 0.0000001); 
 
    x=reshape(x,  Nx, Ny );

   
   %D_r=reshape(D_r,  Nx, Ny );


   %ISNR_Bayes=10*log10(norm(fcolor(20:Nx-20,20:Ny-20)-g__2(20:Nx-20,20:Ny-20),'fro')^2/(norm(fcolor(20:Nx-20,20:Ny-20,1)-x_r(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,2)-x_g(20:Nx-20,20:Ny-20),'fro')^2+norm(fcolor(20:Nx-20,20:Ny-20,3)-x_b(20:Nx-20,20:Ny-20),'fro')^2))

   boundary=decFactor;
    


   ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )

   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1) - z(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )


   SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));
   


  %  Nx=Nx+2*extras;
  %  Ny=Ny+2*extras;
   % Hcubic=fftshift(Hcubic);
 %HF= abs(Hcubic(1:Nx/2,1:Ny/2)).^2+abs(Hcubic(Nx/2+1:Nx,1:Ny/2)).^2+abs(Hcubic(1:Nx/2,Ny/2+1:Ny)).^2+abs(Hcubic(Nx/2+1:Nx,Ny/2+1:Ny)).^2;
    %'alpha denumerator'
   % TrH= sum(sum( HF./(HF+ssigma*alpha*ones(size(HF)))))  ;
   %alpha= Nx*Ny/ ( TrH+alpha*sum(sum(Qp)) + sum(sum( (x-z).^2 ) ) ) 
    
%Hcubic=fftshift(Hcubic);
% Nx=Nx-2*extras;
    %Ny=Ny-2*extras;

end

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
% iter
% if norm(x-x_prev,'fro')<0.0001
%     break;
% end
toc
end
toc
 
%ISNR=10*log10(  ( norm( f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)  -   g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2   )/(  norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-x(boundary:Nx-boundary+1,boundary:Ny-boundary+1),'fro')^2))

%save(sprintf('file%dwithISNR%d',myinput1),'ISNR_Bayes');
 %   PSNR=10*log10(  (Nx-2*(boundary-1)+1)*(Ny-2*(boundary-1)+1)  /  (   norm(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1)-g2(boundary:Nx-boundary+1,boundary:Ny-boundary+1),   'fro'    )^2 ) )
    
 %       SSIMres=ssim(f(boundary:Nx-boundary+1,boundary:Ny-boundary+1,1),gather(x(boundary:Nx-boundary+1,boundary:Ny-boundary+1)));

