function [loglikelihood,nu,c2]=trainPCSingleModel(expNum,x,c2,nu,coord,E,S)

[Nx,Ny]=size(x);
        loglikelihood=-inf;



        for emiter=1:1220
             
            
            ZALL=zeros(size(x));
            
            for exper=1:expNum
                J=(S{1}/2+nu{exper}/2);
                A{exper}  =   gpuArray((nu{exper}+S{1})./(c2{exper}*E{exper}+nu{exper}));
                %  Z{exper}  =  gpuArray( exp(  log(c2{exper})/2+(S{1}/2+nu{exper}/2-1)*(log(A{exper})-log(J)+psi(J))-A{exper}.*(c2{exper}*E{exper}/2+nu{exper}/2) -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2)  ));
              %Z{exper}  =  gpuArray(  PW{exper}.*((gamma((nu{exper}+1)/2)/gamma((nu{exper})/2)) *(c2{exper}/nu{exper})^0.5*(1+c2{exper}*E{exper}/nu{exper} ).^(-(nu{exper}+1)/2 )));%
               % ZALL=ZALL+Z{exper};
                %ZexperZALL=Z{exper}./ZALL;
                %c2{exper}=sum(sum(ZexperZALL))/sum(sum(ZexperZALL.*A{exper}));
                %ZALL=ZALL-Z{exper};
            end

            %prevlikelihood  =   loglikelihood;
            %loglikelihood  =   0;

            improv=0;
            
            for exper=1:expNum
                    J=(S{1}/2+nu{exper}/2);

                   % Z{exper}  =   (Z{exper}+10^(-16)) ./ (ZALL+10^(-16));
                   % PW{exper}=(Z{exper}+0.001)./(1+0.001);

                    c2prev=c2{exper};
                    c2{exper}  =   Nx*Ny  /   sum(sum( A{exper}.*E{exper}   ));
                    improv=improv+abs(c2prev-c2{exper}).^2;
                    nuprev=nu{exper};
                    nu{exper}=bisection_nu(nu{exper},200,A{exper},ones(Nx,Ny),Nx*Ny,S{1});
                    improv=improv+abs(nuprev-nu{exper}).^2;

                    % loglikelihood  =   loglikelihood  +   0.5*sum(sum( -(c2{exper}  +   nu{exper})*A{exper}.*Z{exper}.*E{exper}  +   (nu{exper}-1)*(Z{exper}).*(log(A{exper})-log(J)+psi(J)) )) ;
                    % loglikelihood  =   loglikelihood  +   sum(sum(  Z{exper}   ))*(log(c2{exper})/2 -log(gamma(nu{exper}/2))  +   log(nu{exper}/2)*(nu{exper}/2));
            end
nu{1:expNum}
c2{1:expNum}
             if improv<10^(-8)
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
        
%disp('End')
         