function [nu,c2]=learnPCMixture(expNum,x,c2,nu,Hcubic,HQ,ssigma,S,rw,decFactor,extras)

[Nx,Ny]=size(x);
    Imask=ones(Nx,Ny);
    Imask(1:decFactor:Nx,1:decFactor:Ny)=0;
    Imask= not(Imask);
    hcubatrous=fftshift(real(ifft2( Hcubic )));
    hcubatrous=hcubatrous(extras+1:Nx+extras,extras+1:Ny+extras);
    hcubatrous=hcubatrous.*Imask;
    
    Qp=zeros(Nx,Ny);

    
    
        for emiter=1:2320
             
            ZALL=zeros(Nx,Ny);
            
            
for exper=1:expNum

        J=(S{1}/2+nu{exper}/2);

        
        %E1{exper}=(f-circshift(f,rw{2*exper-1})  ).^2;
      
         E{exper}   =  real(  ifft2( fft2( ( x-circshift(x, rw{2*exper-1})  ).^2  + Qp+circshift(Qp, rw{2*exper-1})) .* conj(HQ)  ));
      %  E   =  real(  ifft2( fft2( ( x_r-circshift(x_r, rw{2*exper-1})  ).^2  +  (x_g-circshift(x_g, rw{2*exper-1})  ).^2 + (x_b-circshift(x_b, rw{2*exper-1})  ).^2 ) .* conj(HQ{exper})  ));

        %E = (real(ifft2( fft2( ( x-circshift(x,rw{2*exper-1})  ).^2  + Qp) .*( HQ)  )));
        A{exper}  =   (nu{exper}+S{1})./(c2{exper}*E{exper}+nu{exper})   ;
        %Z{exper}  =   exp(  J*log(A{exper})-A{exper}.*(c2{exper}*E/2+nu{exper}/2)   );
        %Z{exper}  =   exp(  log(c2{exper})/2+(S{1}/2+nu{exper}/2-1)*(log(A{exper})-log(J)+psi(J))-A{exper}.*(c2{exper}*E/2+nu{exper}/2) -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2)  );
            
        Z{exper}  =    ((gamma((nu{exper}+S{1})/2)/gamma((nu{exper})/2)) *(c2{exper}/nu{exper})^0.5*(1+c2{exper}.*E{exper}/nu{exper} ).^(-(nu{exper}+S{1})/2 ));%


        %c2prev{exper}=c2{exper};

        %if iter >1 
        %c2{exper}=Nx*Ny/sum(sum( (B{exper}/c2prev{exper}).*( x-circshift(x, rw{2*exper-1})  ).^2));
        %end

        %Z{exper}=(nu+c2{exper}*E).^(-J);
        %B{exper}=A{exper};
        ZALL=ZALL+Z{exper};

     end

    for exper=1:expNum
        Z{exper}=Z{exper}./ZALL;
       % PW{exper}=(Z{exper}+5.001)/(1+5.001);
       % pof{exper}=sum(sum(Z{exper}));
    end
      Zen=zeros(Nx,Ny);
    %ZALL=zeros(Nx,Ny);
    for exper=1:expNum
        Zen=Zen-log(Z{exper}).*Z{exper}/(Nx*Ny); 
    end
    
    zenind=find(Zen>-log(1/expNum)/(Nx*Ny)-10^(-6));
    for exper=1:expNum

        Z{exper}(zenind)=0;
       %ZALL=ZALL+Z{exper};
    end
    

    
 
    Qp=zeros(Nx,Ny);
    Qp= Qp+sum(sum(  hcubatrous  .*   hcubatrous))  /    ssigma;
    Qp=Qp.*Imask;
    
    for exper=1:expNum

        Qp=Qp+B{exper}+circshift(ssigma*B{exper},rw{2*exper-1});%real(ifft2(conj(fft2(q1)).*(fft2(B{exper}))));
    end

    Qp=1./Qp;

    
    
          

            improv=0;
            
            for exper=1:expNum
                   J=(S{1}/2+nu{exper}/2);

                   Z{exper}  =   (Z{exper}+10^(-12)) ./ (ZALL+10^(-12));

                   PW{exper}=(Z{exper}+5.001)./(5+0.001);

                   c2prev=c2{exper};
                    
                   c2{exper}  = (  sum(sum(  Z{exper}   )))  /  ( sum(sum(  Z{exper}.*A{exper}  .*E{exper}   )));
                    
                   if c2{exper}>10^5
                       %disp 'Too big lambda'
                       c2{exper}=10^5;
                       %disp 'Sum of Zs:'
                       sum(sum(  Z{exper}   ))
                   end
                   
                   improv=improv+abs(c2prev-c2{exper}).^2;
                   nuprev=nu{exper};
                   nu{exper}=gather(bisection_nu(nu{exper},300,A{exper},Z{exper},Nx*Ny,S{1}));
                   %   nu{exper}=bisection_nu(nu{exper},200,A{exper},ones(Nx,Ny),Nx*Ny,S{1});
                   %   nu{exper}=bisection_nu(nu{exper},200,A{exper},ones(Nx,Ny)+(A{exper}-ones(Nx,Ny)).*Z{exper},Nx*Ny,S{1});
                   improv=improv+abs(nuprev-nu{exper}).^2;
                   % loglikelihood  =   loglikelihood  +   0.5*sum(sum( -(c2{exper}  +   nu{exper})*A{exper}.*Z{exper}.*E{exper}  +   (nu{exper}-1)*(Z{exper}).*(log(A{exper})-log(J)+psi(J)) )) ;
                   % loglikelihood  =   loglikelihood  +   sum(sum(  Z{exper}   ))*(log(c2{exper})/2 -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2));
            end

            (c2{1:expNum})
            (nu{1:expNum})
            improv

            if improv/expNum<10^(-4)
                improv
                disp('End')
                break;
            end

            % if abs(prevlikelihood  - loglikelihood)<0.01
%                  disp('Exit due to loglikelihood convergence ')
%                  prevlikelihood
%                  loglikelihood
%                  emiter
                 %abs(prevlikelihood  - loglikelihood);

             %    break;
         %    end

        end
        
disp('End')
emiter
         