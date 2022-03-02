clear all
close all

% molecola: m*q''=F
% p'=F
% q' = p'

N_molecole=3;
q=zeros (2,N_molecole);
p=zeros (2,N_molecole);
q(:,1)=1/sqrt(2)*[-1,-1];
p(:,1)=[0,0];
q(:,2)=1/sqrt(2)*[1,1];
p(:,2)=[0,0];
q(:,3)=1/sqrt(2)*[-1,1];
p(:,3)=[0,0];

i=1;
for a=1:length(q)
    z0(i:i+1,1)=q(:,a);
    i=i+2;
end
for a=1:length(q)
    z0(i:i+1,1)=p(:,a);
    i=i+2;
end

t0=0;
tf=100;
h=0.001;

fun=@(t,z) system_force(z);

% [yy,nevals,tt]= euler_esplicito (fun, t0, tf, z0, h);
% [yy1,nstep,nrej,cost,Hr,Ha,STIMA,tt1]= RKembedded (fun,t0,tf,z0,@RK_2and3,TOL);
[tt1, yy1,nevals_rk]= RKclassico (fun, t0, tf, h, z0, Verner6);
% [yy1,nevals, tt1] =Storme_Verlet (fun, t0, tf, z0, h);
% plot(yy(:,1),yy(:,2),"*",yy(:,3),yy(:,4),"*",yy(:,5),yy(:,6),"*")
hold on
grid on
yy1=yy1';
plot(yy1(:,1),yy1(:,2),"*",yy1(:,3),yy1(:,4),"*",yy1(:,5),yy1(:,6),"*")
quiver(yy1(:,1),yy1(:,2),yy1(:,7),yy1(:,8))
quiver(yy1(:,3),yy1(:,4),yy1(:,9),yy1(:,10))
quiver(yy1(:,5),yy1(:,6),yy1(:,11),yy1(:,12))

N_molecole=3;
for a=1:length(tt1)
    K_media(a)=0;
    for i=N_molecole*2+1:2:N_molecole*4
       K_media(a)=K_media(a)+0.5*(yy1(a,i).^2+yy1(a,i+1).^2);
    end
end

for a=1:length(tt1)
    b=0;
    for i=1:2:N_molecole*2
        b=b+1;
        q_tt(1,b)=yy1(a,i);
        q_tt(2,b)=yy1(a,i+1);
    end
    U(a)=U_potenziale(q_tt);
end
figure()

plot(tt1,K_media,"*")
hold on
grid on
plot(tt1,U,"*")

figure
semilogy(tt1,abs((K_media+U)-U(1)), "*")
grid on
