function [loglikelihood,nu,c2]=trainPCMixtureOneL(expNum,c2,nu,E,S)
S{1}=1;
[Nx,Ny]=size(E{1});
        loglikelihood=-inf;
        %expNum=expNum;

% 
         for exper=1:expNum
%             E{exper}  =   real(  ifft2( fft2(  (x-circshift(x, coord(exper,:)  )).^2   )  .* conj(HQ{1})   ));
                PW{exper}=ones(Nx,Ny)/expNum;
         end

        for emiter=1:200
             
            ZALL=zeros(Nx,Ny);
            
            for exper=1:expNum
                J=(1/2+nu{exper}/2);
                A{exper}  =   gpuArray((nu{exper}+1)./(c2{exper}*E{exper}+nu{exper}));
                %  Z{exper}  =  gpuArray( exp(  log(c2{exper})/2+(S{1}/2+nu{exper}/2-1)*(log(A{exper})-log(J)+psi(J))-A{exper}.*(c2{exper}*E{exper}/2+nu{exper}/2) -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2)  ));
                Z{exper}  = gpuArray(  ((gamma((nu{exper}+1)/2)/gamma((nu{exper})/2)) *(c2{exper}/nu{exper})^0.5*(1+c2{exper}*E{exper}/nu{exper} ).^(-(nu{exper}+1)/2 ))   );%
                %Z{exper}  = gpuArray(   ((1+c2{exper}*E{exper}/nu{exper} ).^(-(nu{exper}+S{1})/2 )))  ;%
                %Gaussian formula
                %Z{exper}  =  gpuArray(  PW{exper}.*( c2{exper}^0.5*exp(-c2{exper}*E{exper}/2 )) ); %
                ZALL=ZALL+Z{exper};
                %ZexperZALL=Z{exper}./ZALL;
                %c2{exper}=sum(sum(ZexperZALL))/sum(sum(ZexperZALL.*A{exper}));
                %ZALL=ZALL-Z{exper};
            end

            %prevlikelihood  =   loglikelihood;
            %loglikelihood  =   0;

            improv=0;
            lambda=0;

            nuprev=nu{exper};

            mysumA=  (psi((nuprev+S{1})/2) -  log((nuprev+S{1})/2) +1);
                for exper=1:expNum
       Z{exper}=Z{exper}./ZALL;
    %   PW{exper}=(Z{exper}+.1)/(1+.1);
       % pof{exper}=sum(sum(Z{exper}));
    end
    Zen=zeros(Nx,Ny);
    %ZALL=zeros(Nx,Ny);
    for exper=1:expNum
        Zen=Zen-log(Z{exper}).*Z{exper}/(Nx*Ny); 
    end
    
    zenind=find(Zen>-log(1/expNum)/(Nx*Ny)-10^(-7));
    for exper=1:expNum

        Z{exper}(zenind)=0;
       %ZALL=ZALL+Z{exper};
    end
    [z1,z2]=size(zenind);
            for exper=1:expNum
                  %  J=(S{1}/2+nu{exper}/2);
                   % Z{exper}  =   (Z{exper}+10^(-16)) ./ (ZALL+10^(-16));
                 % PW{exper}=(Z{exper}+.001)./(1+.001);
                    c2prev=c2{exper};
                   %c2{exper}  =  gather(   sum(sum(  Z{exper}.*A{exper}  .*E{exper}   )));
                   lambda  = lambda + gather(   sum(sum(  Z{exper}.*A{exper}  .*E{exper}   )));
                    s1=0;
                    %while  s1<Nx*Ny/100
                %        EZ=c2{exper}*Z{exper}.*E{exper}/nu{exper};
                  %      [s1, s2]=size( find(EZ(:)>1) );
                  %      c2{exper}=c2{exper}*(1-s1/(Nx*Ny));
                    %end
            %c2{exper}=max(nu{1}/max(max(E.*Z{exper})),1000);
        
%                     if c2{exper}>10^8
%                         %disp 'Too big lambda'
%                         c2{exper}=10^8;
%                         %disp 'Sum of Zs:'
%                         sum(sum(  Z{exper}   ))
%                     end

                mysumA  =   mysumA+sum(sum(  Z{exper}.*log(A{exper}) - Z{exper}.*A{exper} ))/(Nx*Ny-z1*z2) ;

                  improv=improv+abs(c2prev-c2{exper}).^2;
                  %  nuprev=nu{exper};
                  %  nu{exper}=gather(bisection_nu(nu{exper},300,A{exper},Z{exper},Nx*Ny,S{1}));
                   %   nu{exper}=bisection_nu(nu{exper},200,A{exper},ones(Nx,Ny),Nx*Ny,S{1});
                      %nu{exper}=bisection_nu(nu{exper},200,A{exper},ones(Nx,Ny)+(A{exper}-ones(Nx,Ny)).*Z{exper},Nx*Ny,S{1});
                  %  improv=improv+abs(nuprev-nu{exper}).^2;
                    % loglikelihood  =   loglikelihood  +   0.5*sum(sum( -(c2{exper}  +   nu{exper})*A{exper}.*Z{exper}.*E{exper}  +   (nu{exper}-1)*(Z{exper}).*(log(A{exper})-log(J)+psi(J)) )) ;
                    % loglikelihood  =   loglikelihood  +   sum(sum(  Z{exper}   ))*(log(c2{exper})/2 -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2));
            end

        %  nuall=gather(bisectioAlln_nu(nu{exper},300,mysumA,Nx*Ny,S{1}))
            
            S{1}*gather( Nx*Ny  / lambda )

            for exper=1:expNum
                c2{exper}  =   S{1}*gather( (Nx*Ny-z1*z2)  / lambda );
      %        nu{exper}=nuall;
            end

            %(c2{1:expNum})

            %(nu{1:expNum})
            %improv

%             if improv/expNum<10^(-4)
%                 improv
%                 disp('End')
%                 break;
%             end

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
         