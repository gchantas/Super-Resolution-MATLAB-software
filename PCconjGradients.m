function x=PCconjGradients(x0,fncHandle,b,iw,rw, preCon, maxiter,rtol)

decFactor=iw{2};
Hcubic=iw{3};
extras=iw{4};
Nx=iw{5};
Ny=iw{6};

x=x0;
r=b-fncHandle(x0,iw,rw);
 

% 
% 
% bVec=myimresize(reshape(r,[Nx Ny]),Nx,Ny,Hcubic,decFactor,extras);
% 
% 
% x1=real(ifft2(fft2(bVec).*(psInv)));
% x1_=gpuArray2(zeros(Nx,Ny));
% x1_(1:decFactor:Nx,1:decFactor:Ny)=x1;
% z=decFactor^2*real(ifft2(conj(Hcubic).*fft2(x1_)));
%z=z(:);

 z=preCon.*r;

p = z;

for iter=1:maxiter
    

    Ap=fncHandle(p,iw,rw);    

    a=sum(z.*r)/sum( p.*Ap);

    x=x+a*p;

    rnormprev=sum(r.*z);   

    %if mod(iter,50)~=0
        r=r-a*Ap;
    %else
%       r=b-fncHandle(x,iw,rw);        
%    end
%         disp 'rnorm'
     %rnorm=sum(r.^2)

     
     a*norm(p,'fro')^2
    
     if a*norm(p,'fro')^2<rtol
         disp 'Exit due to convergence: rnorm is lesser than rtol'
         disp 'Iterations made:  '
         iter
         return;
     end

% 
% bVec=myimresize(reshape(r,[Nx Ny]),Nx,Ny,Hcubic,decFactor,extras);
% 
% 
% x1=real(ifft2(fft2(bVec).*(psInv)));
% x1_=gpuArray2(zeros(Nx,Ny));
% x1_(1:decFactor:Nx,1:decFactor:Ny)=x1;
% z=decFactor^2*real(ifft2(conj(Hcubic).*fft2(x1_)));
% z=z(:);
    z=preCon.*r;

    b=sum(z.*r)/rnormprev;

    
    p=z+b*p;
end

disp 'Exit due to reaching the maximum number of iterations '


