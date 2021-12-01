clear all
close all

% y˙= Ay+g(t) y(0) = [2,3/2]' , t ∈ [0, T]
% Soluzione esatta
% syms x(t) y(t)
% B = [2*sin(t); 2*(cos(t)-sin(t))];
% Y = [x;y];
% C = Y(0) == y0;
% odes = diff(Y) == A*Y + B;
% [y1Sol(t),y2Sol(t)] = dsolve(odes, C);
% y1Sol(t) = simplify(y1Sol(t));
% y2Sol(t) = simplify(y2Sol(t));

% Prima A e g(t)
y1_es= @(t)(5*exp(-t))/4 + (3*exp(-3*t))/4 + sin(t);
y2_es= @(t)(5*exp(-t))/4 - (3*exp(-3*t))/4 + cos(t);

% Secondi A e g(t)
% y1_es= @(t)(5987*exp(-t))/1998 + (1004009*exp(-1000*t))/1998001998 - (4*124501124501^(1/2)*cos(t + atan(249751/249250)))/1000001
% y2_es= @(t)(5987*exp(-t))/1998 - (501000491*exp(-1000*t))/999000999 - (2*495013495013^(1/2)*cos(t + atan(497503/497498)))/1000001
t0=0;
tf=10;
y0=[2,3/2]';

A=[-2 1; 1 -2];
%A=[-2 1; 998 -999];

g=@(t) [2*sin(t); 2*(cos(t)-sin(t))];
%g=@(t) [2*sin(t); 999*(cos(t)-sin(t))];

h=0.8;
tol=10e-8;
maxiter=1000;

fun=@(t,y) A*[y(1);y(2)]+g(t);
jac=@(t,y) A;
% sol_esatta=@(t) exp(lambda*t)*(y0 - lambda^2/(lambda^2 + 1)) - (lambda*(sin(t) - lambda*cos(t)))/(lambda^2 + 1);

for d=1
    [yy_Ei,nevals_Ei, tt_Ei]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,0);
    [yy_Ee,nevals_Ee, tt_Ee]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,1);
    [yy_tr,nevals_tr, tt_tr]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,1/2);
    [yy_pm,nevals_pm, tt_pm]=punto_medio(fun,jac,t0, tf , y0, h,tol,maxiter);

    [yy_rke,nstep, nrej, neval, Hr, Ha,STIMA,tt_rke]= RKembedded (fun,t0,tf,y0,@RK_2and3,1e-5);
    [tt_rk,yy_rk, nevals] = RKclassico (fun, t0, tf, h, y0, Verner6);
    %Soluzione esatta t,y
    figure()
    subplot(1,2,1)
    plot(tt_Ei,yy_Ei(1,:), "*-", tt_Ee,yy_Ee(1,:), "*-",tt_pm, yy_pm(1,:), "*-", ...
         tt_tr, yy_tr(1,:), "*-", tt_rke, yy_rke(1,:), "*-", tt_rk, yy_rk(1,:), "*-");
    grid on; hold on; title(['Soluzione 1 con h ', num2str(h)]); xlabel("t"); ylabel("y");
    fplot(y1_es,[t0,tf]);
    legend("Euler implicito", "Euler esplicito", "Punto medio",  ...
        "trapezi","Rk embedded 2-3","Rk esplicito", "Soluzione esatta")
    subplot(1,2,2)
    plot(tt_Ei,yy_Ei(2,:), "*-", tt_Ee,yy_Ee(2,:), "*-",tt_pm, yy_pm(2,:), "*-", ...
         tt_tr, yy_tr(2,:), "*-", tt_rke, yy_rke(2,:), "*-", tt_rk, yy_rk(2,:), "*-");
    grid on; hold on; title(['Soluzione 2 con h ', num2str(h)]); xlabel("t"); ylabel("y");
    fplot(y2_es,[t0,tf]);
    legend("Euler implicito", "Euler esplicito", "Punto medio",  ...
        "trapezi","Rk embedded 2-3","Rk esplicito", "Soluzione esatta")
end

% Studiamo l'h_min per la stabilità
hh=[4,2,1,0.5,0.25,0.125,0.675,0.3,0.15,0.07,0.035];
i=1;
for h=hh
    [yy_Ei,nevals_Ei, tt_Ei]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,0);
    err_Ei(i)= max(norm((yy_Ei(1,:))'-y1_es(tt_Ei),"inf"), norm((yy_Ei(2,:))'-y2_es(tt_Ei),"inf"));
    nevals_Ei_err(i)=nevals_Ei*2/3;
    [yy_pm,nevals_pm, tt_pm]=punto_medio(fun,jac,t0, tf , y0, h,tol,maxiter);
    err_pm(i)= max(norm((yy_pm(1,:))'-y1_es(tt_pm),"inf"), norm((yy_pm(2,:))'-y2_es(tt_pm),"inf"));
    nevals_pm_err(i)=nevals_pm;
    [yy_Ee,nevals_Ee, tt_Ee]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,1);
    err_Ee(i)= max(norm((yy_Ee(1,:))'-y1_es(tt_Ee),"inf"), norm((yy_Ee(2,:))'-y2_es(tt_Ee),"inf"));
    nevals_Ee_err(i)=nevals_Ee*2/3;
    [yy_tr,nevals_tr, tt_tr]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,0.5);
    err_tr(i)= max(norm((yy_tr(1,:))'-y1_es(tt_tr),"inf"), norm((yy_tr(2,:))'-y2_es(tt_tr),"inf"));
    nevals_tr_err(i)=nevals_tr;
    [tt_rk,yy_rk, nevals] = RKclassico (fun, t0, tf, h, y0, Heun3);
    err_rk(i)= max(norm((yy_rk(1,:))'-y1_es(tt_rk),"inf"), norm((yy_rk(2,:))'-y2_es(tt_rk),"inf"));
    nevals_rk_err(i)=nevals;

    i=i+1;
end
figure()
title(['Errore vs h' ])
subplot(2,2,1, 'XScale', 'log', 'YScale', 'log')
grid on; hold on; title(['Eulero' ]); xlabel("h"); ylabel("err");
loglog(hh,err_Ei, "*", hh,err_Ee, "*");
legend("Euler implicito", "Euler esplicito");

subplot(2,2,2, 'XScale', 'log', 'YScale', 'log')
grid on; hold on; title(['Trapezi' ]); xlabel("h"); ylabel("err");
loglog(hh, err_tr, '*');
legend("trapezi");

subplot(2,2,3, 'XScale', 'log', 'YScale', 'log')
grid on; hold on; title(['Punto medio' ]); xlabel("h"); ylabel("err");
loglog(hh, err_pm, '*');
legend("Punto medio");

subplot(2,2,4, 'XScale', 'log', 'YScale', 'log')
grid on; hold on; title(['rk ordine6' ]); xlabel("h"); ylabel("err");
loglog(hh, err_pm, '*');
legend("Rk6");

% loglog(hh,err_Ei, "*", hh,err_Ee, "*",hh, err_pm, '*',hh, err_tr, '*',hh, err_rk, '*');
% grid on; hold on; title(['Errore vs h' ]); xlabel("h"); ylabel("err");
% legend("Euler implicito", "Euler esplicito", "Punto medio","trapezi", "rk V6")

% dt al variare del tempo
figure()
semilogy(Ha(:,1), Ha(:,2), "*")
hold on
semilogy(Hr(:,1), Hr(:,2), "*")
legend(['Accettati Tool ', num2str(tol)], ['Rifiutati Tool ', num2str(tol)]);
grid on; title('t vs h'); xlabel("t"); ylabel("h");

% nevals vs errore 
figure()
loglog(err_Ei,nevals_Ei_err, "*",err_Ee,nevals_Ee_err, "*", ...
    err_pm,nevals_pm_err, "*",err_tr,nevals_tr_err, "*",err_rk,nevals_rk_err, "*");
legend("Euler implicito", "Euler esplicito", "Punto medio","trapezi", "rk V6")
grid on; title('t vs h'); xlabel("t"); ylabel("h");title(['Diagramma di efficienza' ])


