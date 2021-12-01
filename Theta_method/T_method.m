function [yy,nevals, tt]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,theta,rpar)

if (size(y0,2)~=1)
     y0=y0';
end

nevals=0;
N = ceil ((tf-t0)/h)+1; 
tt=zeros(N,1);
tt(1)=t0;
yy=zeros(size(y0,1),N);
yy(:,1)=y0;

k = zeros(1,N-1);
dim= size(y0,1);

G = @(Yk,Yn,h,t_aux) Yk - h*( (1-theta)*fun(t_aux+h,Yk) + theta*fun(t_aux,Yn)) - Yn;
Gp = @(Yk,h,t_aux) eye(dim) - h*(1-theta)*jac(t_aux,Yk);

for n=1:N-1
    Yn = yy(:,n);
    Yk = Yn;
    err = 1;

    if (tt(n)+h>tf)
        h=tf-tt(n);
    end

    while (err>tol)&&(k(n)<maxiter)
        Ykp1 = Yk - Gp(Yk,h,tt(n))\G(Yk,Yn,h,tt(n));
        err = norm(Ykp1-Yk);
        Yk = Ykp1;
        k(n) = k(n)+1;
        nevals=nevals+3; %1 jac e 2 fun
    end

    yy(:,n+1) = Ykp1;
     tt(n+1)=tt(n)+h;    
end
