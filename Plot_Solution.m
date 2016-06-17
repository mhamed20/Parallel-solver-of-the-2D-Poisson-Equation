clear all

 m=15;
 n=15;

 load x.dat;
 load y.dat;
%  load f.dat
 load dx.dat;
 load ddx.dat;

 % Define grid
for i=1:m
    for j=1:n      
         X(i,j)=x(i+(j-1)*m,1);
         Y(i,j)=y(i+(j-1)*m,1);
%          F(i,j)=f(i+(j-1)*m,1);
      end
end          
   for i=1:m
    for j=1:n          
            Du(i,j)=dx(i+(j-1)*m,1);
            DDu(i,j)=ddx(i+(j-1)*m,1);
    end
   end
   
% 
%    for i=2:m-1
%     for j=2:n-1   
%      x_xi=(X(i+1,j)-X(i-1,j))/2;
%      y_xi=(Y(i+1,j)-Y(i-1,j))/2;
%      u_xi=(u[i+1+ m*(j)]-u[i-1 + m*(j)])/2;
%      x_etha=(X(i,j+1)-)/2;
%      y_etha=(y[i + m*(j+1)]-y[i + m*(j-1)])/2;
%      u_etha=(u[i + m*(j+1)]-u[i + m*(j-1)])/2;
%      J=(x_xi*y_etha-x_etha*y_xi);
%      D0(i,j)=(1/J)*(u_xi*y_etha-u_etha*y_xi);     
%     end
%    end
%  
% mesh(dx,x,y)
% mesh(ddx,x,y)

    