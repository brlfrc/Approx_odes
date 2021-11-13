% clear all
% close all

% Pendolo: m*L*(theta)''=-m*g*sin(theta(t))
% E(t)=1/2*m*(L*theta')^2-m*L*g*cos(theta(t))
% p = theta'

g= 9.8;
L= 1;
m= 1;

for N=[25,50,100,200,400]
    t0=0;
    tf=2;
    h=tf/N;
    y0=[pi/4; 0]; %[theta, p]
    
    f=@(t,y) [y(2); -g/L*sin(y(1))];

    [yy,nevals,tt]= euler_esplicito (f, t0, tf, y0, h);
    
    figure(1)
    plot(tt,yy(:,1),"*")
    grid on
    hold on
end

figure(1)
title("#theta di un pendolo")
legend("N 25", "N 50", "N 100", "N 200", "N 400")

T=1/2*m*(L*yy(:,2)).^2;
U=-m*L*g*cos(yy(:,1));
figure(2)
title("Energie")
plot(tt,T+U,"*",tt,T,"*",tt,U,"*")
grid on
legend("E totale", "E cinetica", "E potenziale")

E_esatta=U(1);
figure(3)
plot(tt, abs(E_esatta-(T+U)), "*")
title("Errore Energia meccanica, con N=400")
grid on
hold on
xlabel("t")
ylabel("E_meccanica")

%Non ho fatto vedere, ma al variare di n aumenta la pendenza della curva abs(E_esatta-(T+U)) 
