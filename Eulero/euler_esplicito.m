function[y,nevals, tt]= euler_esplicito (fun, t0, tf, y0, h)

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

for j=1:N
    
    yy(j+1,:)=yy(j,:) +h*fun(t_aux,yy(j,:))'; %FUNZIONI passate come vettori in colonna
    nevals= nevals+1;

    if (t_aux+h>tf)
        h=tf-t_aux;
    end
    t_aux=t_aux+h;
    tt(j+1)=t_aux;

    
end


y=yy;
end

%Così funziona ma non mi piace moltissimo, perché con un intervallo esatto
%mette un doppio punto finale che fa un po' schifo, ma alemeno ora
%funziona. In generale avrei preferito un metodo che prende in input il
%numero di intervalli e non la larghezza di questi, ma va beh.