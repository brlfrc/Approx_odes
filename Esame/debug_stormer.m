clear all
close all

% molecola: m*q''=F
% p'=F
% q' = p'

N_molecole=100;
q=10*rand(2,N_molecole);
p=10*rand(2,N_molecole);


i=1;
for a=1:length(q)
    z0(i:i+1,1)=q(:,a);
    i=i+2;
end
for a=1:length(q)
    z0(i:i+1,1)=p(:,a);
    i=i+2;
end

t0=1;
tf=1.1;
h=0.0001;
fun=@(t,z) system_force(z);

[yy,nevals, tt] = Storme_Verlet (fun, t0, tf, z0, h);


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

plot(tt,K_media,"*")
hold on
plot(tt,U,"*")

figure
semilogy(tt,K_media+U, "*")

% for i=1:2:N_molecole*2
%     hold on
%     plot(yy(:,i),yy(:,i+1),"*");
% end
% legend();
% plot(yy(:,1),yy(:,2),"*",yy(:,3),yy(:,4),"*",yy(:,5),yy(:,6),"*")
% hold on
% toc;
