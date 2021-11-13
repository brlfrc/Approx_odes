function [y,nstep, nrech, nevals, H_r,STIMA,tt]= RKembedded (f,t0,tf,y0,tableau,TOL)

[c,A,b,B,p,P] = feval(tableau);

h0 = max(TOL^(1/(p+1)),0.25);

n=length(y0);

tt=zeros(1, 1);
tt(1)=t0;
t_aux=t0;


yy=zeros(n, 1);
yy(:,1)=y0;
y_aux=y0;
y_aux0=y0;
y_aux1=y0-5*TOL;

s= size(A,1);
K=zeros(n, s); % [ K1 | K2| K3...] vettori colonna

cost=0;
m=0;
nstep=0;
nrech=0;
j=0;

while t_aux<tf
    h=h0;
    count=0;
    while max(abs(y_aux1-y_aux0))>TOL
        est=max(abs(y_aux1-y_aux0));
        count=count+1;
        h=h*0.75*max(0.4,(TOL/(est+0.01))^(1/(min(P,p)+1)));
        temp=0;
    
        for i=1:s
            for j=1:i
                temp = temp + A(i,j)*K(:,j);
            end
            
            K(:,i) = f(t_aux + c(i)*h, y_aux + temp*h);
            cost=cost+1;
            temp=0;
        end
        temp1=0;
        for i=1:s
                temp = temp +  b(i)*K(:,i);
                temp1= temp1 + B(i)*K(:,i);
        end
        
        y_aux0 = y_aux+h*temp;
        y_aux1 = y_aux+h*temp1;

        nrech= nrech+1;
        H_r(1,j)=t_aux;
        H_r(2,j)=h;
        j=j+1;

        if(count>1000)
            break;
        end
   end
    y_aux= y_aux+h*temp;
    yy(:,m+1)=y_aux;

    if (t_aux+h>tf)
        h=tf-t_aux;
    end
    t_aux=t_aux+h;
    tt(m+1)=t_aux;
    m=m+1;

    nstep=nstep+1;
    STIMA(1,nstep)=t_aux;
    STIMA(2,nstep)=max(abs(y_aux1-y_aux));
    y_aux1=y_aux-5*TOL;
end

nevals=cost;
y=yy;
nrech=nrech-nstep;

end
