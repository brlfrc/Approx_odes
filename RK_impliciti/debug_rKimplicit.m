clear all
close all

y0=[0.001];
lambda=3;
fun=@(t,y) [lambda*y*(1-y)];
jac=@(t,y) [lambda*(1-y)-lambda*y];
t0=0;
tf=5;
h=0.1;
%RK_implicit(fun,jac,t0,tf,y0,tableau,maxiter,h,TOL)
[yy,nevals,iter,nfail,tt]=RK_implicit(fun,jac,t0,tf,y0,@GL4,100,h,10^-5);

sol_esatta=@(t)(y0*exp(lambda*t))/(y0*exp(lambda*t) - y0 + 1);

plot(tt,yy,"*")%,tt,yy(2,:),"*")
hold on
fplot(sol_esatta,[t0,tf])
legend("y1", "y2")
grid on;