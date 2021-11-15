clear all
close all

% Pendolo: m*L*(theta)''=-m*g*sin(theta(t))
% E(t)=1/2*m*(L*theta')^2-m*L*g*cos(theta(t))
% p = theta'

g= 9.8;
L= 1;
m= 1;

h = 10e-3;
TOL =10e-5;

t0=0;
tf=2;
y0=[pi/4; 0]; %[theta, p]

f=@(t,y) [y(2); -g/L*sin(y(1))];

[yy_EU,nevals_EU,tt_EU]= euler_esplicito (f, t0, tf, y0, h);
[tt_RK,yy_RK, nevals_RK] = RKclassico (f, t0, tf, h, y0, Fehlberg4);
[yy_ERK,nstep, nrech, nevals_ERK, H_r,STIMA,tt_ERK]= RKembedded (f,t0,tf,y0,@Fehlberg_4and5,TOL);

% Digramma soluzioni
figure(1)
title("soluzione")
plot(tt_EU,yy_EU(:,1),"*", tt_RK,yy_RK(1,:),"*", tt_ERK,yy_ERK(1,:),"*")
xlabel("t"); ylabel("#theta")
grid on
legend("Eulero", "RK classico", "RK embedded");

%Diagramma sulle Energie

E_EU= 1/2*m*(L*yy_EU(:,2)).^2 -m*L*g*cos(yy_EU(:,1));
E_RK= 1/2*m*(L*yy_RK(2,:)).^2 -m*L*g*cos(yy_RK(1,:));
E_ERK= 1/2*m*(L*yy_ERK(2,:)).^2 -m*L*g*cos(yy_ERK(1,:));
E_esatta=E_RK(1);

figure()
semilogy(tt_EU, abs(E_esatta-E_EU), "*", tt_RK, abs(E_esatta-E_RK), "*",tt_ERK, abs(E_esatta-E_ERK), "*")
title("Errore Energia meccanica")
legend("Eulero", "RK classico", "RK embedded");
grid on
xlabel("t")
ylabel("E_meccanica")
