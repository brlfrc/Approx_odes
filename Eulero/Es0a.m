clear all
close all

for h=[0.25,0.1, 0.05]

    t0=0;
    tf=5;
    y0=1;
    [yy,nevals, tt]= euler_esplicito (@(t,y) -5*y, t0, tf, y0, h);
    hold on;
    plot(tt,yy,"*");

end

fplot(@(x) exp(-5*x), [0,5]);
legend("h=0.25","h=0.1","h=0.05", "sol esatta");