clear all
close all

% y'=cos^2(ay), y in R^n e t in [0,T]
alpha=40;
g=@(t,y) cos(alpha*y)*cos(alpha*y);

y0=-0;
c=-tan(alpha*y0)/alpha;
f_esatta=@(t) atan(alpha*(t-c))/alpha;

t0=0;
tf=0.4;
h=10e-3;
TOL=10e-7;

[yy_E,nevals_E,tt_E]= euler_esplicito (g, t0, tf, y0, h);
[yy,nstep, nrech, nevals, H_r,H_a,STIMA,tt]= RKembedded (g,t0,tf,y0,@EulerHeun  ,TOL);

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
n=size(STIMA, 1);
STIMA=STIMA';
semilogy(STIMA(1, 1:n-1), abs(STIMA(1, 1:n-1)-STIMA(1,2:n)), "*", tt_E, h*ones(length(tt_E)), "*")
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
