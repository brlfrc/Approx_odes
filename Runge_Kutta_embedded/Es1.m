clear all
close all

% y'=y^2, y in R^n e t in [0,T]
% RISOLTO IN MODO COMPLETO

t0=0;
y0=2;
tf=0.45;
h=10e-3;
TOL=10e-7;

g=@(t,y) y*y;
f_esatta=@(t) y0./(1-t*y0);

[yy_E,nevals_E,tt_E]= euler_esplicito (g, t0, tf, y0, h);
[yy,nstep, nrech, nevals, H_r, H_a,STIMA,tt]= RKembedded (g,t0,tf,y0,@RK_2and3,TOL);

%Soluzione esatta t,y
figure()
title(['Soluzione con toll = ', num2str(TOL)])
xlabel("t"); ylabel("h");

plot(tt,yy, "*", tt_E,yy_E, "*");
grid on
hold on
fplot(f_esatta, [t0,tf])
legend("ERK", "Euler", "sol esatta")

% dt al variare del tempo
figure()
semilogy(H_a(:,1), H_a(:,2), "*")
hold on
semilogy(H_r(:,1), H_r(:,2), "*")
legend(['Accettati Tool ', num2str(TOL)], ['Rifiutati Tool ', num2str(TOL)]);
grid on; title('t vs h'); xlabel("t"); ylabel("h");

% Stima dell'errore locale data dal metodo
figure()

semilogy(STIMA(:, 1), STIMA(:,2), "*")
legend(['toll ', num2str(TOL)])
grid on; title('t vs STIMA ERR loc'); xlabel("t"); ylabel("ERR loc");

% Errore ERR locale al variare del tempo (sono solo quelli accettati)
figure()
i=1;
n=length(tt);
for t=tt(1:n-1)
    [yy,nstep, nrech, nevals, H_r, H_a,STIMA,tt_aa]= RKembedded (g,t,tf,f_esatta(t),@RK_2and3,TOL);
    err_local(i)=abs(f_esatta(H_a(1,1)+H_a(1,2))-yy(1,2));
    i=i+1;
end
semilogy(tt(1:n-1),err_local, "*")
legend(['toll ', num2str(TOL)])
grid on
title('t vs ERR loc')
xlabel("t");
ylabel("ERR loc");
