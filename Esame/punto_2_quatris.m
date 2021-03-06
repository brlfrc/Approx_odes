clear all
close all

N_molecole=100;
for i=0:9
    q(:,(i+1)+9*i:(i+1)+9*i+9)=[i+zeros(1,10);
        0:9];
end
p=zeros(2,N_molecole);
z0=convert(q,p);
t0=1;
tf=2;
h=[0.1,0.01,0.001,0.0001];
fun=@(t,z) system_force(z);
fun1=@(t,z) system_force_modified(z);

for j=1:length(h)
    j
    [yy,nevals,tt, count1] = Storme_Verlet(fun, t0, tf, z0, h(j));
    diff_norm(j)=en(tt,yy,N_molecole);

    [yy,nevals,tt, count2] = Storme_Verlet(fun1, t0, tf, z0, h(j));
    diff_modified(j)=en(tt,yy,N_molecole);
end

figure
loglog(h(1:4),diff_norm, "*",h(1:4),diff_modified, "*")
xlabel("h")
ylabel("E_M")
grid on
legend("StormerVerlet classico", "StormerVerlet modificato")

function z0= convert(q,p)
    i=1;
    for a=1:length(q)
        z0(i:i+1,1)=q(:,a);
        i=i+2;
    end
    for a=1:length(q)
        z0(i:i+1,1)=p(:,a);
        i=i+2;
    end
end

function diff= en(tt,yy,N_molecole)
    for a=1:length(tt)
        K_media(a)=0;
            for i=N_molecole*2+1:2:N_molecole*4
               K_media(a)=K_media(a)+0.5*(yy(a,i).^2+yy(a,i+1).^2);
            end
    end
    
    for a=1:length(tt)
        b=0;
        for i=1:2:N_molecole*2
            b=b+1;
            q_tt(1,b)=yy(a,i);
            q_tt(2,b)=yy(a,i+1);
        end
        U(a)=U_potenziale(q_tt);
    end

    diff=max(abs(K_media+U-K_media(1)-U(1)));
end