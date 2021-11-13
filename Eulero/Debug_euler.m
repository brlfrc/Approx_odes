clear all
close all

%Debug Euler_esplicito con h non multiplo dell'intervallo

t0= 0;
tf= 1;
h= 0.25;
y0=0;

[y,nevals, tt]= euler_esplicito (@(t,y) 1, t0, tf, y0, h);

plot(tt,y,"*");
hold on
grid on
fplot(@(x) x, [0,1])