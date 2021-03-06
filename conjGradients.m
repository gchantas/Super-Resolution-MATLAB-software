function x=conjGradients(x0,fncHandle,b,iw,rw,  maxiter,rtol)

x=x0;
r=b-fncHandle(x0,iw,rw);

p = r;
rnormStart=sum(r.^2);
rnorm=sum(r.^2);
for iter=1:maxiter
    
%iter

%         if mod(iter,50)==0
%             r=b-fncHandle(x,iw,rw);
% 
%             p=r;
%         end
    Ap=fncHandle(p,iw,rw);    

    a=sum(r.^2)/sum( p.*Ap);

    
    
       %if rnorm<rtol*rnormStart
       if mod(iter,20)==0
            if a * norm(p,'fro')^2/norm(x,'fro')^2   <     rtol
              %  disp 'Exit due to convergence: rnorm is lesser than rtol'
              %  disp 'Iterations made:  '
              %  iter
                return;
            end
       end
    
    

        
    x=x+a*p;

    

%           if mod(iter,50)==0
%              r=b-fncHandle(x,iw,rw);
%                            p=r;
%                            continue
%           else
           
              r=r-a*Ap;
       %   end
    
    rnormprev=rnorm;  

    
        
    rnorm=sum(r.^2);

 

    beta=rnorm/rnormprev;

    
    p=r+beta*p;
    

end

disp 'Exit due to reaching the maximum number of iterations '


