function c2=bisection_lamda_gt(c2,uplimit,A,N,M,nprev,c)

NN=N*M; %NN=N*N;

up=uplimit;
down=c2;


first = grad_gtn(down,A,nprev,c,N, M);
sec   = grad_gtn(up,A,nprev,c,N, M);

if first*sec>0
    c2=-1;
    return;
end;
    
while abs(first)>1
    
first =  grad_gtn(down,A,nprev,c,N, M);

mid=(up+down)/2;



third =  grad_gtn(mid,A,nprev,c,N, M);

if third*first>0
    down=mid;
else
    up=mid;
end

end

c2=mid;