function vec=tageDiagBlockMat(A,Nx,Ny,decFac,iter)

if decFac==1
   
    vec=A;
else

    vec=gpuArray2(zeros(Nx,Ny));
        
        vec(1:Nx/2,1:Ny/2)=tageDiagBlockMat(A{1,1}{1,1},Nx/2,Ny/2,decFac/4,iter); 
        vec(Nx/2+1:Nx,1:Ny/2)=tageDiagBlockMat(A{1,1}{2,2},Nx/2,Ny/2,decFac/4,iter); 
        vec(Nx/2+1:Nx,Ny/2+1:Ny)=tageDiagBlockMat(A{2,2}{1,1},Nx/2,Ny/2,decFac/4,iter); 
        vec(1:Nx/2,Ny/2+1:Ny)=tageDiagBlockMat(A{2,2}{2,2},Nx/2,Ny/2,decFac/4,iter);

end



