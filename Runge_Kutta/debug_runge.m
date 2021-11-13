clear all
close all

fun= @(t,y) [y(1) ; 1];
t0=0;
tf=1;
y0=[1;1];
h=0.0013;

[tt,yy, nevals] = RKclassico (fun, t0, tf, h, y0, Verner6);

plot(tt,yy(1,:), "*",tt,yy(2,:), "*")
grid on
hold on
fplot(@(x) exp(x), [t0,tf])
fplot(@(x) x+1, [t0,tf])

