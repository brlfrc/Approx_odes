clear all
close all

a=1;

lambda=-10;
    
t0=0;
tf=5;
y0=1;

for k=[4,5,6]
    TOL=10^(-k);
    [yy,nstep, nrech, nevals, H_r,STIMA,tt]= RKembedded (@(t,y) lambda*y,t0,tf,y0,@EulerHeun,TOL);
    n=size(STIMA, 2);
    
    figure(1)
    loglog(STIMA(1, 1:n-1), abs(STIMA(1, 1:n-1)-STIMA(1,2:n)), "*")
   
    hold on
end

% h vs dt
title(['lambda = ', num2str(lambda)])
grid on
xlabel("t"); ylabel("h");
legend("toll 10e-4", "toll 10e-5", "toll 10e-6")


figure()
title(['lambda = ', num2str(lambda), ' Toll= 10e-',  num2str(k)])

hold on;
grid on

plot(tt,yy,"*");
sol_esatta= @(x) exp(lambda*x);
fplot(sol_esatta, [0,5]);
legend("soluzione approssimata","soluzione esatta");
