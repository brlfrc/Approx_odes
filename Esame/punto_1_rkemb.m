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
tf=200;
h=0.001;
TOL=1e-2;
TOL1=1e-6;
fun=@(t,z) system_force(z);

[yy,nstep1,nrej1,cost,Hr,Ha,STIMA,tt]= RKembedded (fun,t0,tf,z0,@Fehlberg_4and5,TOL);
yy=yy';
[yy1,nstep, nrej, nevals, Hr1, Ha1,STIMA1,tt1]= RKembedded (fun,t0,tf,z0,@RK_2and3,TOL1);
yy1=yy1';

figure()
semilogy(Ha1(1:50:end,1), Ha1(1:50:end,2), "*")
hold on
semilogy(Hr1(:,1), Hr1(:,2), "*")
legend(['Accettati Tool ', num2str(TOL1)], ['Rifiutati Tool ', num2str(TOL1)]);
grid on; xlabel("t"); ylabel("h");

figure()
semilogy(Ha(1:100:end,1), Ha(1:100:end,2), "*")
hold on
semilogy(Hr(:,1), Hr(:,2), "*")
legend(['Accettati Tool ', num2str(TOL)], ['Rifiutati Tool ', num2str(TOL)]);
grid on; xlabel("t"); ylabel("h");

figure()
semilogy(STIMA(1:100:end, 1), STIMA(1:100:end,2), "*")
legend(['Stima errore locale con toll ', num2str(TOL)])
grid on; xlabel("t"); ylabel("ERR loc");

figure()
semilogy(STIMA1(1:100:end, 1), STIMA1(1:100:end,2), "*")
legend(['Stima errore locale con toll ', num2str(TOL1)])
grid on; xlabel("t"); ylabel("ERR loc");

% hold on
% grid on
% plot(yy1(:,1),yy1(:,2),"b*",yy1(:,3),yy1(:,4),"y*",yy1(:,5),yy1(:,6),"r*")
% quiver(yy1(:,1),yy1(:,2),yy1(:,7),yy1(:,8))
% quiver(yy1(:,3),yy1(:,4),yy1(:,9),yy1(:,10))
% quiver(yy1(:,5),yy1(:,6),yy1(:,11),yy1(:,12))
% xlabel("x")
% ylabel("y")

% N_molecole=3;
% for a=1:length(tt1)
%     K_media(a)=0;
%     for i=N_molecole*2+1:2:N_molecole*4
%        K_media(a)=K_media(a)+0.5*(yy1(a,i).^2+yy1(a,i+1).^2);
%     end
% end
% 
% for a=1:length(tt1)
%     b=0;
%     for i=1:2:N_molecole*2
%         b=b+1;
%         q_tt(1,b)=yy1(a,i);
%         q_tt(2,b)=yy1(a,i+1);
%     end
%     U(a)=U_potenziale(q_tt);
% end
% figure()
% 
% plot(tt1,K_media,"*")
% hold on
% grid on
% plot(tt1,U,"*")
% xlabel("t")
% ylabel("E")
% legend("Cinetica", "Potenziale")
% 
% figure
% semilogy(tt1,abs((K_media+U)-U(1)), "*")
% grid on
% xlabel("t")
% ylabel("\DeltaE(t)")
% 
% figure
% plot(tt1,abs((K_media+U)-U(1))/abs(U(1)), "*")
% grid on
% xlabel("t")
% ylabel("\DeltaE_%(t)")


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