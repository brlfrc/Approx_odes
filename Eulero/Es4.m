clear all
close all

% y'=-y^(-1/2) in R^n e t in [0,T]

t0=0;
tf_1=1/2;
tf_2=3/2;
h=0.01;
y0=1;
g=@(t,y) -y^(-1/2);

[yy_1,~,tt_1]= euler_esplicito (g, t0, tf_1, y0, h);
[yy_2,nevals,tt_2]= euler_esplicito (g, t0, tf_2, y0, h);
[t_ode45,y_ode45] = ode45(g, [0,3/2], y0);

plot(tt_1,yy_1, "*-", tt_2,yy_2, "*",t_ode45,y_ode45);
legend("Tfinale=1/2", "finale=3/2", "ode45")

% Osservo che non gunziona più. Perchè? Cade Lipschtiz?