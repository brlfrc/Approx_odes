% Faccio leggermente diverso rispetto a quanto proposto
% Non riesco a capire come passare in argomento la funzione sarà da
% sistemare



function [tt,yy, nevals] = RKclassico (fun, t0, tf, h, y0, tableau)

nstep= ceil ((tf-t0)/h);
n=length(y0);

tt=zeros(nstep+1, 1);
tt(1)=t0;
t_aux=t0;

yy=zeros(n, nstep+1);
yy(:,1)=y0;
y_aux=y0;

AA= tableau;
n_AA=size(AA,1);
c= AA(1:n_AA-1,1);
b= AA(n_AA,2:n_AA)';
A= AA(1:n_AA-1,2:n_AA);

s= size(A,1);
K=zeros(n, s); % [ K1 | K2| K3...] vettori colonna

cost=0;

for m=1:nstep
    
    temp=0;

    for i=1:s
        for j=1:i
            temp = temp + A(i,j)*K(:,j);
        end
        
        K(:,i) = fun(t_aux + c(i)*h, y_aux + temp*h);
        cost=cost+1;
        temp=0;
    end
    
    for i=1:s
            temp = temp +  b(i)*K(:,i); %si può rendere matriciale?
    end
    y_aux = y_aux+h*temp;

    yy(:,m+1)=y_aux;

    if (t_aux+h>tf)
        h=tf-t_aux;
    end
    t_aux=t_aux+h;
    tt(m+1)=t_aux;

end

nevals=cost;

end
