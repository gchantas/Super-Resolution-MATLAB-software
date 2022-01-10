function nu=bisection_nu(nu,uplimit,A,nuprev,NN,S)



up=uplimit;
down=nu;

%mysumA=sum(sum( Z.*log(A./Z) ))-sum(sum(A))+ sum(sum( Z*(psi((nu2+S)/2)-log((nu+S)/2)) - 2*(1-Z)/(nu+S)  ));
mysumA=sum(sum( log(A) ))/NN-sum(sum(A))/NN+ sum(sum( (psi((nuprev+S)/2)-log((nuprev+S)/2))   ));
first = mysumA  - ((psi(down/2)-1-log(down/2)));
sec   = mysumA  - ((psi(up/2)-1-log(up/2)));

if first*sec>0
    nu=-1;
    return;
end;
    
while abs(first)>1
    
first =  mysumA  - NN*((psi(down/2)-1-log(down/2)));

mid=(up+down)/2;



third =  mysumA   - NN*((psi(mid/2)-1-log(mid/2)));

if third*first>0
    down=mid;
else
    up=mid;
end

end

nu=mid;