clear all
close all

i=1;
t0=0;
tf=5;

y0=0;
lambda=-2;

hh_c=1./[2,4,8,16,32];

sol_esatta=@(t) exp(lambda*t)*(y0 - lambda^2/(lambda^2 + 1)) - (lambda*(sin(t) - lambda*cos(t)))/(lambda^2 + 1);
fun=@(t,y) lambda*(y-cos(t));
jac=@(t,y) lambda;
for h=hh_c
    [yy_pm,nevals_tr, tt_pm]=punto_medio(fun,jac,t0, tf , y0, h,10^(-6),1000);
    err_pm_c(i)= norm(yy_pm'-sol_esatta(tt_pm),"inf");
    [yy_tr,nevals_tr, tt_tr]=T_method(fun,jac,t0, tf , y0, h,10^(-6),1000,0.5);
    err_tr_c(i)= norm(yy_tr'-sol_esatta(tt_tr),"inf");
    i=i+1;
end
for n=2:size(err_pm_c,2)
    ord=err_tr_c(n)/err_tr_c(n-1)
end

figure()
loglog(hh_c, err_pm_c, '*');
grid on; hold on; title(['Errore vs h. A lambda= ', num2str(lambda) ]); xlabel("h"); ylabel("err");
legend("Punto medio")