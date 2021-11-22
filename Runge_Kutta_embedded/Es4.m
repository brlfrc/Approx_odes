clear all
close all

y0=[3/2; 3];
f=@(t,y) [1+y(1)^2*y(2)-4*y(1); 3*y(1)-y(1)^2*y(2)];

t0=0;
tf=20;

TOL_1=10e-5;
TOL_2=10e-8;


% Digramma delle fasi ala variare del punto iniziale
figure()
for k= [1,2,3,4,5]
    y0=[k/2; k+2];
    [yy,tt,nstep,nrej,cost,est,Ha,Hr]= RKembedded (f,t0,tf,y0,@RK_2and3,TOL_1);
    plot(yy(1,:),yy(2,:), "*");
    hold on
end
grid on; title('Piano delle fasi'); xlabel("y1"); ylabel("y2");

[yy_1,nstep, nrech, nevals, H_r_1,H_a_1,STIMA_1,tt_1]= RKembedded (f,t0,tf,y0,@RK_2and3,TOL_1);
[yy_2,nstep, nrech, nevals, H_r_2,H_a_2,STIMA_2,tt_2]= RKembedded (f,t0,tf,y0,@RK_2and3,TOL_2);

%Soluzione a diverse tollerenza
figure()
subplot(1,2,1)
plot(tt_1,yy_1(1,:), "*",tt_2,yy_2(1,:), "*");
legend(['Tool ', num2str(TOL_1)], ['Tool ', num2str(TOL_2)])
grid on; title('Soluzione con diverse toll y1'); xlabel("t"); ylabel("y1");

subplot(1,2,2)
plot(tt_1,yy_1(2,:), "*",tt_2,yy_2(2,:), "*");
legend(['Tool ', num2str(TOL_1)], ['Tool ', num2str(TOL_2)]);
grid on; title('Soluzione con diverse toll y2'); xlabel("t"); ylabel("y2");

% dt al variare del tempo
figure()
semilogy(H_a_1(:,1), H_a_1(:,2), "*", H_a_2(:,1), H_a_2(:,2), "*")
hold on
semilogy(H_r_1(:,1), H_r_1(:,2), "*", H_r_2(:,1), H_r_2(:,2), "*")
legend(['Accettati Tool ', num2str(TOL_1)], ['Accettati Tool ', num2str(TOL_2)], ...
       ['Rifiutati Tool ', num2str(TOL_1)], ['Rifiutati Tool ', num2str(TOL_2)]);
grid on; title('t vs h'); xlabel("t"); ylabel("h");

% Stima dell'errore locale data dal metodo
figure()

semilogy(STIMA_1(:, 1), STIMA_1(:,2), "*",STIMA_2(:, 1), STIMA_2(:,2), "*")
legend(['toll ', num2str(TOL_1)], ['toll ', num2str(TOL_2)])
grid on; title('t vs STIMA ERR loc'); xlabel("t"); ylabel("ERR loc");