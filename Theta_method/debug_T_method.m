clear all
close all

fun= @(t,y) [y(1); 1];
t0=0;
tf=1;
y0=[1;1];
h=0.001;
tol= 10e-8;
maxiter=100;
jac=@(t,y) [1,0;0,0];
theta=0.5;

[yy,nevals, tt]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,theta);

plot(tt,yy(1,:), "*",tt,yy(2,:), "*")
grid on
hold on
fplot(@(x) exp(x), [t0,tf])
fplot(@(x) x+1, [t0,tf])

