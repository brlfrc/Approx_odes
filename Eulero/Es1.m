clear all
close all

a=1;

for lambda= i*1 %[-1,-10,-100,1,10]
    
    t0=0;
    tf=5;
    y0=1;
    h=0.01;
    [yy,nevals, tt]= euler_esplicito (@(t,y) lambda*y, t0, tf, y0, h);

    figure()
    title(["lambda = ", num2str(lambda)])
    
    hold on;
    grid on

    plot(tt,yy,"*");
    sol_esatta= @(x) exp(lambda*x);
    fplot(sol_esatta, [0,5]);
    legend("soluzione approssimata","soluzione esatta");

    err_conv(a)= max(abs(sol_esatta(tt)-yy));
    a=a+1;
end


figure()
title(["Errore del metodo"])
semilogy([-1,-10,-100,1,10],err_conv,"*");
grid on
