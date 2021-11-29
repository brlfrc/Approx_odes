clear all
close all

%y'=10y(1-y)

t0=0;
y0=0.0001;
tf=2;
h=10e-3;

g=@(t,y) 10*y*(1-y);
jac=@(t,y) 10-20*y;


tau=0.1;
n=1;
tbar=0.95;
sol_esatta=@(t) 1./(1+exp(-(t-tbar)/tau));

[yy,nevals,tt]= punto_medio(g,jac,t0, tf , y0, h,10e-8,1000);
[yy_E0,nevals_E0,tt_E0]= T_method(g,jac,t0, tf , y0, h,10e-8,1000,0);
[yy_E1,nevals_E1,tt_E1]= T_method(g,jac,t0, tf , y0, h,10e-8,1000,1);

%Soluzione esatta t,y
figure()
plot(tt_E0,yy_E0, "*", tt_E1,yy_E1, "*",tt,yy, "*");
grid on;xlabel("t"); ylabel("h");title("Logistica con eulero implicito e esplicito")
hold on
fplot(sol_esatta,[t0,tf]);
legend("Euler explicit", "Euler implicit", "Punto medio", "Sol esatta")
