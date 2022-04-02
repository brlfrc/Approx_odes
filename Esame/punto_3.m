clear all
close all
% Fissimo il numero di particelle 100

% q=grid_hexagonal(10000,1.1);
% plot(q(1,:),q(2,:),".", "Markersize",25)
% grid on
% xlabel("x")
% ylabel("y")


% i=0;
% l=[0.9:0.01:1.2];
% for L=l
%     i=i+1;
%     U(i)=U_potenziale(grid_hexagonal(10000,L));
% end
% 
% figure
% plot(l,U,"*")
% grid on;
% xlabel("l") 
% ylabel("U")
% 0.98 1.107 1.2
% 1.11
N_molecole=10000;
q=grid_hexagonal(N_molecole,1.107);
p=zeros(2,N_molecole);

z0=convert(q,p);
t0=0;
tf=3;
h=0.001;
fun=@(t,z) system_force(z);
fun1=@(t,z) system_force_modified(z);

[yy,nevals,tt] = Storme_Verlet(fun1, t0, tf, z0, h);


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

figure
plot(tt,K_media,"*")
hold on
plot(tt,U,"*")
grid on
xlabel("t")
ylabel("E")
legend("Cinematica", "Potenziale")

figure
semilogy(tt,abs(K_media+U-K_media(1)-U(1)), "*")
grid on
xlabel("t")
ylabel("\DeltaE(t)")

figure
for i=1:2:N_molecole*2
    hold on
    plot(yy(:,i),yy(:,i+1),"*");
end
xlabel("x")
ylabel("y")
grid on

figure
plot(tt,abs((K_media+U)-U(1))/abs(U(1)), "*")
grid on
xlabel("t")
ylabel("\DeltaE_%(t)")

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



