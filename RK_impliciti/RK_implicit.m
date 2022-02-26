function [yy,nevals,iter,nfail,tt]=RK_implicit(fun,jac,t0,tf,y0,tableau,maxiter,h,TOL)
if (size(y0,2)~=1)
     y0=y0';
end

[c,A,b] = feval(tableau);

iter=0;
nfail=0;
nevals=0;

N = ceil ((tf-t0)/h)+1; 
tt=zeros(N,1);
tt(1)=t0;
yy=zeros(size(y0,1),N);
yy(:,1)=y0;

s=length(c);
dim= size(y0,1);

F= @(z,h,tn,yn,c,fun) -z + h*kron(A,eye(dim))*F_Z(z,h,tn,yn,c,fun);
dF = @(z,h,Jac) -eye(dim*s) + h*kron(A,Jac);  % per ora problema monodimensionale poi devo cambiare qui

for n=1:N-1

    Jac=jac(tt(n),yy(:,n));
    z=zeros(length(c),1);
    Z0 = Z_vector(z,dim,s);

    Z = newton(F, dF, Z0, TOL, maxiter,yy(:,n),h,tt(n),c,fun, Jac);
    Z=reshape(Z,[dim,s]);
    Z=Z';
    yy(:,n+1) = yy(:,n)+(b'*(A\Z))';
    tt(n+1)=tt(n)+h;
    if (tt(n)+h>tf)
        h=tf-tt(n);
    end
end

% iter=sum(k);
% nfail=sum(k)-N-1;

end

function Z = Z_vector(z,dim,s)
    Z=zeros(dim,s);
    for i=1:s
        Z(:,i)=z(i);
    end
    Z=Z(:);
end

function f_Z = F_Z(z,h,tn,yn,c,fun) %passo z matrice
    dim=length(yn);
    s=length(c);
    Z=reshape(z,[dim,s]);
    f_Z= zeros(dim,s);
    for j=1:s
            f_Z(:,j)=fun(tn+c(j)*h, yn+Z(:,j));
    end
    f_Z=f_Z(:);
end