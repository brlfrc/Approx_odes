clear all
close all

% y'=y^2, y in R^n e t in [0,T]

t0=0;
y0=0.03;
tf=1/y0;
h=0.5;

g=@(t,y) y*y;
f_esatta=@(t) y0/(1-t*y0);

[yy,nevals,tt]= euler_esplicito (g, t0, tf, y0, h);

tspan = [0 tf];
[t_ode45,y_ode45] = ode45(g, tspan, y0);
[t_ode23,y_ode23] = ode23s(g, tspan, y0);

plot(tt,yy, "*", t_ode45,y_ode45, "*",t_ode23,y_ode23, "*");
grid on
hold on
fplot(f_esatta, [t0,tf])
legend("Soluzione approx", "sol ode45", "sol ode23", "sol esatta")

%Manca convergenza