clear all
close all

% y'=y^2, y in R^n e t in [0,T]

t0=0;
y0=2;
tf=1/y0-0.2;
h=10e-3;
TOL=10e-5;

g=@(t,y) y*y;
f_esatta=@(t) y0./(1-t*y0);

[yy_E,nevals_E,tt_E]= euler_esplicito (g, t0, tf, y0, h);
[yy,nstep, nrech, nevals, H_r,STIMA,tt]= RKembedded (g,t0,tf,y0,@EulerHeun,TOL);

%Soluzione esatta t,y
figure()
title(['Soluzione con toll = ', num2str(TOL)])
xlabel("t"); ylabel("h");

plot(tt,yy, "*", tt_E,yy_E, "*");
grid on
hold on
fplot(f_esatta, [t0,tf])
legend("RK", "Euler", "sol esatta")

% dt al variare del tempo (sono solo quelli accettati)
figure()
n=size(STIMA, 2);
hh=ceil(2:n/50:n);
xx=STIMA(1, hh);
zz=abs(STIMA(1,hh-1)-STIMA(1,hh));

semilogy(xx, zz, "*", tt_E, h*ones(size(tt_E,1)), "*")
legend("toll 10e-5", "Euler")
grid on
title('t vs h')
xlabel("t");
ylabel("h");

% Errore ERR al variare del tempo (sono solo quelli accettati)
figure()
semilogy(tt,abs(yy-f_esatta(tt)), "*", tt_E,abs(yy_E-f_esatta(tt_E)), "*")
legend("toll 10e-5", "Euler")
grid on
title('t vs ERR')
xlabel("t");
ylabel("ERR");
