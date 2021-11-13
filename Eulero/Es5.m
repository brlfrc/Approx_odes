clear all
close all

% y'=cos^2(ay)

a=40;
t0=0;
y0=0;%-pi/(2*a);
tf=2;
h=0.01;

g=@(t,y) (cos(a*y))^2;
f_esatta=@(t) atan(a*(t))/a;  %è sbagliata, manca c che renda verla prima condizione (t-c) dato iniziale !=0

[yy,nevals,tt]= euler_esplicito (g, t0, tf, y0, h);

tspan = [t0 tf];
[t_ode45,y_ode45] = ode45(g, tspan, y0);
[t_ode23,y_ode23] = ode23s(g, tspan, y0);

plot(tt,yy, "*", t_ode45,y_ode45, "*",t_ode23,y_ode23, "*");
grid on
hold on
fplot(f_esatta, [t0,tf])
legend("Soluzione approx", "sol ode45", "sol ode23", "sol esatta")

%Manca convergenza