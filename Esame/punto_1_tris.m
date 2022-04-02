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
h=[0.1,0.05,0.01,0.005];
TOL=1e-4;
fun=@(t,z) system_force(z);


for j=1:length(h)
    j
    [tt, yy,nevals]= RKclassico (fun, t0, tf, h(j), z0, Verner6);
    nevals_Ver(j)=nevals;
    [yy1,nevals,tt1]= euler_esplicito (fun, t0, tf, z0, h(j));
    nevals_eu(j)=nevals;
    [tt, yy,nevals]= RKclassico (fun, t0, tf, h(j), z0, Heun3);
    nevals_Heu(j)=nevals;

end

figure
loglog(h,nevals_Ver, "*",h,nevals_eu, "*",h,nevals_Heu, "*")
xlabel("h")
ylabel("nevals")
grid on



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