function [y,nstep, nrech, nevals, H_r,STIMA,tt]= RKembedded (f,t0,tf,y0,tableau,TOL)

[~,~,~,~,p,P] = feval(tableau);

h0 = max(TOL^(1/(p+1)),0.0001);

tt=zeros(1, 1);
tt(1)=t0;
t_aux=t0;

n=length(y0);
yy=zeros(n, 1);
yy(:,1)=y0;
y_aux=y0;

cost=0;
m=0;
nstep=0;
nrech=0;

while t_aux<tf
    h=h0;
    count=0;
    est=10*TOL;

    while 1
        h=h*0.75*max(0.4, (TOL/(est+0.01))^(1/(min(P,p)+1))); %min(2,(TOL/(est+0.01))^(1/(min(P,p)+1)));
        
        [y_aux0, y_aux1,cost_onestep,y_step]= RK_onestep (y_aux, t_aux,h, tableau, n, f);

        cost= cost + cost_onestep;
        count=count+1;
        nrech= nrech+1;
        H_r(1,nrech)=t_aux;
        H_r(2,nrech)=h;

        est=max(abs(y_aux1-y_aux0));
        
        if(count>100000 || est<TOL)
            break;
        end
    end

    if (t_aux+h>tf)
        h=tf-t_aux;
    end

    y_aux= y_aux+h*y_step;
    yy(:,m+1)=y_aux;

    t_aux=t_aux+h;
    tt(m+1)=t_aux;
    m=m+1;

    nstep=nstep+1;
    STIMA(1,nstep)=t_aux;
    STIMA(2,nstep)=max(abs(y_aux1-y_aux));
end

nevals=cost;
y=yy;
nrech=nrech-nstep;

end

function [y_aux0, y_aux1, cost,y_step]= RK_onestep(y_aux, t_aux, h, tableau, n, f)

[c,A,b,B] = feval(tableau);

s= size(A,1);
K=zeros(n, s);
cost=0;
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
        temp1 = temp1 + B(i)*K(:,i);
end
    
y_aux0 = y_aux+h*temp;
y_aux1 = y_aux+h*temp1;

y_step=temp;
end
