clear all
close all

N_molecole=3;
q(:,1)=1/sqrt(2)*[-1,-1];
p(:,1)=[0,0];
q(:,2)=1/sqrt(2)*[1,1];
p(:,2)=[0,0];
q(:,3)=1/sqrt(2)*[-1,1];
p(:,3)=[0,0];

z0=convert(q,p);
t0=0;
tf=20;
h=[0.1,0.01,0.001,0.0001, 0.00001, 0.000001];
TOL=1e-4;
fun=@(t,z) system_force(z);


for j=1:length(h)
    if j<5
        [tt, yy,nevals_rk]= RKclassico (fun, t0, tf, h(j), z0, Verner6);
        yy=yy';
        diff_v6(j)=en(tt,yy,N_molecole);
        
        [tt, yy,nevals_rk]= RKclassico (fun, t0, tf, h(j), z0, Heun3);
        yy=yy';
        diff_H3(j)=en(tt,yy,N_molecole);
    end

    [yy,nevals,tt]= euler_esplicito (fun, t0, tf, z0, h(j));
    diff_eu(j)=en(tt,yy,N_molecole);

end

figure
loglog(h(1:4),diff_v6, "*",h(1:4),diff_H3, "*",h,diff_eu, "*")
xlabel("h")
ylabel("E_M")
grid on
legend("Verner 6", "Heun 3", "Eulero esplicito")



function z0=convert(q,p)
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