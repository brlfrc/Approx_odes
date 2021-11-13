clear all
close all

% Parametri del problema

A=[-2 1; 1 -2];
fun= @(t,y) A* [y(1); y(2)] + [2*sin(t); 2*(cos(t)-sin(t))];
t0=0;
tf=10;
y0=[2;3];
h=1;
for n=0:12
    hh(n+1)=h;
    h=h/2;
end


% Soluzione Esatta
syms x(t) y(t)
B = [2*sin(t); 2*(cos(t)-sin(t))];
Y = [x;y];
C = Y(0) == y0;
odes = diff(Y) == A*Y + B;
[y1Sol(t),y2Sol(t)] = dsolve(odes, C);
y1Sol(t) = simplify(y1Sol(t));
y2Sol(t) = simplify(y2Sol(t));

% Risolviamo il sistema per diversi h, salviamo le informazioni sull'errore
% e sul numero di valutazioni
i=0;
for h=hh
    i=i+1;
    [tt,yy, nevals_RK4(i)] = RKclassico (fun, t0, tf, h, y0, Heun3); %ordine 4
    sol_esatta =double([y1Sol(tt), y2Sol(tt)]');
    err_RK4 (i)= max(max(abs (yy- sol_esatta)));  %Errore nel prendere la norma
    [tt,yy, nevals_RK2(i)] = RKclassico (fun, t0, tf, h, y0, dop5); %ordine 3
    err_RK2 (i)= max(max(abs (yy- sol_esatta)));
    [yy,nevals_euler(i), tt]= euler_esplicito (fun, t0, tf, y0, h);
    err_euler (i)= max(max(abs (yy'- sol_esatta)));
%     [tt,yy] = ode45(fun, [t0 tf], y0);
%     err_ode45 (i)= norm (yy- [y1Sol(tt), y2Sol(tt)]', inf);
end

% Diagramma di convergenza h versus err in scala log log o semilogy
% RISALE MA CREDO SIA OK
figure()
loglog(hh, err_RK4, "*", hh, err_RK2, "*", hh, err_euler, "*")%, hh, err_ode45, "*")
legend("RK_3", "RK_5", "Euler")%, "ode45")
grid on
title("Errore vs h dei metodi")
xlabel("h"); ylabel("err");

% Diagrammi di efficienza costo versus errore 
figure()
loglog(err_RK4, nevals_RK4, "*", err_RK2, nevals_RK2, "*", err_euler, nevals_euler, "*")
legend("RK_3", "RK_5", "Euler")
grid on
title("Efficienza dei metodi");
xlabel("err"); ylabel("nevals");

figure()
title("Diagramma fasi")
plot(double(y1Sol(tt)),double(y2Sol(tt)), "*")
grid on

% Variazione errore locale
h=0.25;
i=0;
for t0 =[1, 2, 3, 4, 5]
    i=i+1;
    y0= double([y1Sol(t0), y2Sol(t0)]');
    [tt,yy, nevals_RK4(i)] = RKclassico (fun, t0, tf, h, y0, Fehlberg4);
    sol_esatta =double([y1Sol(t0+h), y2Sol(t0+h)]');
    err_RK4_loc (i)= norm (yy(2)- sol_esatta(2), inf);
    [tt,yy, nevals_RK2(i)] = RKclassico (fun, t0, tf, h, y0, Runge3);
    err_RK2_loc (i)= norm (yy(2,:)- sol_esatta(2,:), inf);
    [yy,nevals_euler(i), tt]= euler_esplicito (fun, t0, tf, y0, h);
    err_euler_loc (i)= norm (yy(:,2)- sol_esatta(2,:), inf);
end

figure()
title("Erorelocale a diversi t0")
xlabel("t0"); ylabel("err loc");
loglog([1, 2, 3, 4, 5], err_RK4_loc, "*", [1, 2, 3, 4, 5], err_RK2_loc, "*", [1, 2, 3, 4, 5], err_euler_loc, "*")%, hh, err_ode45, "*")
legend("RK_3", "RK_6", "Euler")
grid on