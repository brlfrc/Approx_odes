clear all
close all

for h=[0.1,0.05, 0.01]
    t0=0;
    tf=0.4;
    y0=1.1;
    [yy,nevals, tt]= euler_esplicito (@(t,y) 50*(y-1).^2*(y-5), t0, tf, y0, h);
    hold on;
    plot(tt,yy,"*");
end
grid on
legend("h=0.1","h=0.05","h=0.01");

% Soluzione esatta molto complessa. In generale abbiamo anche che per
% piccoli intervalli il risultato oscilla molto e non esce