clear all
close all


f= @(t,y) [y(1) ; 1];
t0=0;
tf=1;
y0=[1;1];
TOL= 10^-7;
[y,nstep, nrej, nevals, Hr, Ha,STIMA,tt]= RKembedded (f,t0,tf,y0,@RK_2and3,TOL);

plot(tt,y(1,:), "*",tt,y(2,:), "*")
grid on
hold on
fplot(@(x) exp(x), [t0,tf])
fplot(@(x) x+1, [t0,tf])