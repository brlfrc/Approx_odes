function [yy,nevals, tt] = Storme_Verlet (fun, t0, tf, y0, h)
% Storme_Verlet Implementazione One-step formulation del metodo Stoerme
% Verlet. Applicata al caso molecole
%   fun: funzione che descrive la dinamica (q'=p, p'=F)
%   [t0,tf] tempo
%   z0=[q0,p0] condizioni iniziali
%   h passo
%   yy vettore contenente [q(t), p(t)]
%   nevals numero di valutazioni
%   tt vettore dei tempi

if (size(y0,2)~=1)
     y0=y0';
end

nevals=0;
N = ceil ((tf-t0)/h); % Occhio che qui sto lasciando volontariamente indietro un possibile intervallo (esempio: t0=0, tf=1, h=0.3, da 3 ma gli intervalli sono 4)
tt=zeros(N+1,1);
t_aux=t0;
tt(1)=t0;
yy=zeros(N+1,size(y0,1));
yy(1,:)=y0;

N_part= length(y0)/2;

for j=1:N
    
    yy(j+1,:)=0;
    a=fun(t_aux,yy(j,:));
    v= yy(j,N_part+1:end)+(0.5*h*a(N_part+1:end))';
    yy(j+1,1:N_part)= yy(j,1:N_part)+h*v;
    a=fun(t_aux,yy(j+1,:));
    yy(j+1,N_part+1:end)= v+(h/2*a(N_part+1:end))';
    nevals= nevals+2;

    if (t_aux+h>tf)
        h=tf-t_aux;
    end
    t_aux=t_aux+h;
    tt(j+1)=t_aux;    
end
end