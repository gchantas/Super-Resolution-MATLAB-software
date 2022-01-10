function bVec=subfoldVec(vec,decFac,iter)

if decFac==1
   
    bVec=vec;
else


    
    bVec{1}=subfoldVec(vec,decFac/2,iter+1);
    bVec{2}=subfoldVec(vec,decFac/2,iter+1);    
    
    
end