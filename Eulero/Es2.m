clear all
close all

% y'=Ay+g(t), y in R^n e t in [0,T]

t0=0;
tf=10;
h=0.01;
y0=[2;3/2];
A=[-2 1; 1 -2];
g=@(t,y) A*y'+[2*sin(t); 2*(cos(t)-sin(t))];

[yy,nevals,tt]= euler_esplicito (g, t0, tf, y0, h);

plot(tt,yy(:,1), "*", tt, yy(:,2), "*");
grid on
legend("Soluzione y1", "Soluzione y2")