clear all
close all

% y˙= λ(y − cos(t)) y(0) = y0 , λ < 0.

% Soluzione Esatta
% syms y(t) y0 lambda
% C = y(0) == y0;
% odes = diff(y) == lambda*(y-cos(t));
% y1Sol(t) = dsolve(odes, C);
% y1Sol(t) = simplify(y1Sol(t));
%sol_esatta=exp(lambda*t)*(y0 - lambda^2/(lambda^2 + 1)) - (lambda*(sin(t) - lambda*cos(t)))/(lambda^2 + 1)

t0=0;
tf=0.5;

y0=2;
lambda=-2000;

h=1e-3;
tol=10e-8;
maxiter=100;

fun=@(t,y) lambda*(y-cos(t));
jac=@(t,y) lambda;
sol_esatta=@(t) exp(lambda*t)*(y0 - lambda^2/(lambda^2 + 1)) - (lambda*(sin(t) - lambda*cos(t)))/(lambda^2 + 1);

for d=1
    [yy_Ei,nevals_Ei, tt_Ei]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,0);
    [yy_Ee,nevals_Ee, tt_Ee]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,1);
    [yy_tr,nevals_tr, tt_tr]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,1/2);
    [yy_pm,nevals_pm, tt_pm]=punto_medio(fun,jac,t0, tf , y0, h,tol,maxiter);

    [yy_rke,nstep, nrej, neval, Hr, Ha,STIMA,tt_rke]= RKembedded (fun,t0,tf,y0,@RK_2and3,1e-5);
    [tt_rk,yy_rk, nevals] = RKclassico (fun, t0, tf, h, y0, Verner6);
    %Soluzione esatta t,y
    figure()
    plot(tt_Ei,yy_Ei, "*", tt_Ee,yy_Ee, "*",tt_pm, yy_pm, "*", ...
         tt_tr, yy_tr, "*", tt_rke, yy_rke, "*", tt_rk, yy_rk, "*");
    grid on; hold on; title(['Soluzione a lambda ', num2str(lambda), 'e h ', num2str(h)]); xlabel("t"); ylabel("y");
    fplot(sol_esatta, [t0,tf], "y")
    legend("Euler implicito", "Euler esplicito", "Punto medio",  ...
        "trapezi","Rk embedded 2-3","Rk esplicito", "Sol esatta")
end

% Studiamo l'erroe al variare di h
hh= 10.^-(2:0.1:3);%[1,0.5,0.25,0.125,0.0625,0.0312,0.0156,0.0078];con lambda=200%10.^-(1:7);
i=1;

for h=hh
    [yy_Ei,nevals_Ei, tt_Ei]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,0);
    err_Ei(i)= norm(yy_Ei'-sol_esatta(tt_Ei),"inf");
    [yy_Ee,nevals_Ee, tt_Ee]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,1);
    err_Ee(i)= norm(yy_Ee'-sol_esatta(tt_Ee),"inf");
    [yy_pm,nevals_tr, tt_pm]=punto_medio(fun,jac,t0, tf , y0, h,tol,maxiter);
    err_pm(i)= norm(yy_pm'-sol_esatta(tt_pm),"inf");
    [yy_tr,nevals_pm, tt_tr]=T_method(fun,jac,t0, tf , y0, h,tol,maxiter,0.5);
    err_tr(i)= norm(yy_tr'-sol_esatta(tt_tr),"inf");
    [tt_rk,yy_rk, nevals] = RKclassico (fun, t0, tf, h, y0, Verner6);
    err_rk(i)= norm(yy_rk'-sol_esatta(tt_rk),"inf");

    i=i+1;
end

%Errore contro h
figure()
loglog(hh,err_Ei, "*", hh,err_Ee, "*",hh, err_pm, '*',hh, err_tr, '*',hh, err_rk, '*');
grid on; hold on; title(['Errore vs h. A lambda= ', num2str(lambda) ]); xlabel("h"); ylabel("err");
legend("Euler implicito", "Euler esplicito", "Punto medio","trapezi", "rk V6")
