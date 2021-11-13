clear all
close all

% y'=y^2, y in R^n e t in [0,T]

t0=0;
y0=2;
tf=1/y0-0.1;
h=10e-5;
TOL=10e-5;

g=@(t,y) y*y;
f_esatta=@(t) y0/(1-t*y0);

[yy_E,nevals_E,tt_E]= euler_esplicito (g, t0, tf, y0, h);
[yy,nstep, nrech, nevals, H_r,STIMA,tt]= RKembedded (g,t0,tf,y0,@Fehlberg_4and5,TOL);

figure()
title(['Soluzione con toll = ', num2str(TOL)])
xlabel("t"); ylabel("h");

plot(tt,yy, "*", tt_E,yy_E, "*");
grid on
hold on
fplot(f_esatta, [t0,tf])
legend("RK", "Euler", "sol esatta")


figure()
title('h vs t')
xlabel("t"); ylabel("h");
n=size(STIMA, 2);   %bisogna mettere menop dati
loglog(STIMA(1, 1:n-1), abs(STIMA(1, 1:n-1)-STIMA(1,2:n)), "*", tt_E, h*ones(size(tt_E,1)), "*")
legend("toll 10e-5", "Euler")
grid on
