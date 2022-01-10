function [z,Breturn]= VariationalDenoise(z,y,B,Nx,Ny,coord,expNum,alpha,maxX,maxY)


coord(:,1)=mod(coord(:,1)+Nx/2,Nx)-Nx/2;%mod(k1,Nx);
coord(:,2)=mod(coord(:,2)+Ny/2,Ny)-Ny/2;%mod(l1,Ny);

%maxX=8;%max(abs(mod(coord(:,1)+3*Nx,Nx)))
%maxY=8;%max(abs(mod(coord(:,2)+3*Ny,Ny)))
%maxX=0;
%maxY=0;
nx1=(1:Nx)+maxX;
ny1=(1:Ny)+maxY;
%notnx1=[0:maxX-1, Nx-maxX:Nx-1];
%notny1=[0:maxY-1, Ny-maxY:Ny-1];
Ball=gpuArray(zeros(Nx+2*maxX,Ny+2*maxY));



for k=1:expNum
    Ball(nx1,ny1)=Ball(nx1,ny1)+B{k}(nx1,ny1);
end


for k=1:expNum
    Ball(nx1,ny1) = Ball(nx1,ny1) + B{k}(  nx1+coord(k,1) ,   ny1+coord(k,2) );
  %  Ball(notnx1+1,notny1+1) = Ball(notnx1+1,notny1+1) + B{k}(  mod( notnx1+coord(k,1),   Nx)+1  ,   mod(notny1+coord(k,2),  Ny   )+1   );
end


z1=gpuArray(alpha*y);


for k=1:expNum
    z1(nx1,ny1) =z1(nx1,ny1) + B{k}(nx1,ny1).*z(   nx1-coord(k,1), ny1-coord(k,2) );
   % x1(notnx1+1,notny1+1) = x1(notnx1+1,notny1+1) + B{k}(notnx1+1,notny1+1).*x(   mod( notnx1-coord(k,1), Nx)+1 , mod(notny1-coord(k,2)  ,   Ny) +1 );
end


for k=1:expNum
    z1(nx1,ny1) = z1(nx1,ny1) + B{k}(   nx1+coord(k,1), ny1+coord(k,2)).*z(   nx1+coord(k,1), ny1+coord(k,2));
    %x1(notnx1+1,notny1+1) = x1(notnx1+1,notny1+1) + B{k}(  mod( notnx1+coord(k,1), Nx)  +   1,  mod(notny1+coord(k,2)   , Ny)+1)  .*   x(   mod( notnx1+coord(k,1), Nx)+1 , mod(notny1+coord(k,2)  ,   Ny) +1 );
end


z=z1(nx1,ny1)./(alpha+Ball(nx1,ny1));


Breturn=Ball(nx1,ny1);
% Ball=gpuArray(zeros(Nx,Ny));
% 
% 
% 
% for k=1:expNum
%   Ball=Ball+B{k};
% end
% 
% 
% for k=1:expNum
%   Ball = Ball + circshift( B{k}, -coord(k,:) );
% end
% 
% 
% x1=gpuArray(alpha*y);
% 
% 
% for k=1:expNum
%     x1 = x1+circshift( x,coord(k,:)).*(B{k});
% end
% 
% for k=1:expNum
%     x1 = x1+circshift( x.*(B{k}) ,  -coord(k,:));
% end
% 
% x=x1./   (alpha+Ball);