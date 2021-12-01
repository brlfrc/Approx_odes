function [yy,nevals, tt]=punto_medio(fun,jac,t0, tf , y0, h,tol,maxiter,rpar)

if (size(y0,2)~=1)
     y0=y0';
end
nevals=0;
err=1;
N = ceil ((tf-t0)/h)+1;
tt=zeros(N,1);
tt(1)=t0;

yy=zeros(size(y0,1),N);
yy(:,1)=y0;

dim= size(y0,1);

G = @(z,Yn,h,t_aux) -z + h*(fun(t_aux+h/2,(Yn+z)/2)) + Yn;
Gp = @(z,h,t_aux) -eye(dim) + h/2*jac(t_aux,z);

for n=1:N-1
    Yn = yy(:,n);
    Z0 = Yn;

    if (tt(n)+h>tf)
        h=tf-tt(n);
    end
    
    k=0;
    
   [Z,m] = newton(G, Gp, Z0, tol, maxiter,Yn,h,tt(n));
%     while (err>tol)&&(k<maxiter)
%         Ykp1 = Z - Gp(Z,h,tt(n))\G(Z,Yn,h,tt(n));
%         err = norm(Ykp1-Z);
%         Z = Ykp1;
%         k = k+1;
%         nevals=nevals+2; %1 jac e 1 fun
%     end
    nevals=nevals+m;
    yy(:,n+1) = Z;
    tt(n+1)=tt(n)+h;    
end