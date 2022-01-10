function bVec=makeBlockVec(vec,Nx,Ny,decFac,iter)

if decFac==1
   
    bVec=vec;
else

    if mod(iter,2)==0
        D1{1}= vec(1:Nx/2,1:Ny/2); 
        D1{2}= vec(Nx/2+1:Nx,1:Ny/2); 
        D2{1}= vec(Nx/2+1:Nx,Ny/2+1:Ny); 
        D2{2}= vec(1:Nx/2,Ny/2+1:Ny); 
        Nx=Nx/2;
        Ny=Ny/2;
        bVec{1}=makeBlockVec(D1,Nx,Ny,decFac/2,iter+1);
        bVec{2}=makeBlockVec(D2,Nx,Ny,decFac/2,iter+1);
    else
        bVec{1}=makeBlockVec(vec{1},Nx,Ny,decFac/2,iter+1);
        bVec{2}=makeBlockVec(vec{2},Nx,Ny,decFac/2,iter+1);
    end

end
