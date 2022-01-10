function nu=bisection_nu(nu,uplimit,mysumA,A,Z,NN,S)

nuprev=nu;

up=uplimit;
down=.000001;

%mysumA=sum(sum( Z.*log(A./Z) ))-sum(sum(A))+ sum(sum( Z*(psi((nu2+S)/2)-log((nu+S)/2)) - 2*(1-Z)/(nu+S)  ));
%mysumA=sum(sum( Z.*log(A) - Z.*A ))/sum(sum(Z))   +  psi((nuprev+S)/2) -  log((nuprev+S)/2) +1  ;

%mysumA=sum(sum( log(A) - A ))/NN   +  (psi((nuprev+S)/2)-log((nuprev+S)/2))  ;


first = mysumA  -  (psi(  (down)/2)  -  log(  (down)/2))  ;
sec   = mysumA  -  (psi(  (up)/2)  -  log(  (up)/2))  ;

if first*sec>0
    
    nu=up
    disp('Exit due to negative nu')
    
    return;
end;
    
mid=nu;
 
while abs(first)>0.00000001
    
    first =  mysumA  -((psi(down/2)-log(down/2)));

    mid=(up+down)/2;



    third =  mysumA   - ((psi(mid/2)-log(mid/2)));

    if third*first>0
        down=mid;
    else
        up=mid;
    end

    if third==first || third==0 ||first==0
        nu=mid;
        return;
    end
    
end

nu=mid;