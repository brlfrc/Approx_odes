clear all
close all

% y'=(a-by)y, y in R^n e t in [0,T]

t0=0;
tf=5;
y0=5;
a=2.5;
b=1;
g=@(t,y) (a-b*y)*y;
f_esatta=@(t) a*y0/(b*y0+(a-b*y0)*exp(-a*t));

[yy,nevals,tt]= euler_esplicito (g, t0, tf, y0, h);

tspan = [0 5];
[t_ode45,y_ode45] = ode45(g, tspan, y0);
[t_ode23,y_ode23] = ode23s(g, tspan, y0);

figure()
plot(tt,yy, "*", t_ode45,y_ode45, "*",t_ode23,y_ode23, "*");
grid on
hold on
fplot(f_esatta, [t0,tf])
legend("Soluzione approx", "sol ode45", "sol ode23", "sol esatta")

%Manca convergenza