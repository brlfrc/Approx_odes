clear all
close all

g= 9.8;
L= 1;
m= 1;

for N=[25,50,100,200,400]
    t0=0;
    tf=2;
    h=tf/N;
    y0=[pi/4; 0];
 
    f=@(t,y) [y(2); -g/L*sin(y(1))];

    [yy,nevals,tt]= euler_esplicito (f, t0, tf, y0, h);
    [tt_rk, yy_rk,nevals_rk]= RKclassico (f, t0, tf, h, y0, Verner6);
    
end

figure(1)
title("#theta di un pendolo")
plot(tt,yy(:,1),"*")
hold on
plot(tt_rk,yy_rk(1,:),"*")
grid on
legend("euler", "Rk ordine 6")

T=1/2*m*(L*yy(:,2)).^2;
T_rk=1/2*m*(L*yy_rk(2,:)).^2;
U=-m*L*g*cos(yy(:,1));
U_rk=-m*L*g*cos(yy_rk(1,:));
figure(2)
title("Energie eulero")
plot(tt,T+U,"*",tt,T,"*",tt,U,"*")
grid on
legend("E totale", "E cinetica", "E potenziale")

figure(3)
title("Energie rk")
plot(tt_rk,T_rk+U_rk,"*",tt_rk,T_rk,"*",tt_rk,U_rk,"*")
grid on
legend("E totale", "E cinetica", "E potenziale")

E_esatta=U(1);
figure(4)
plot(tt, abs(E_esatta-(T+U)), "*", tt_rk, abs(E_esatta-(T_rk+U_rk)), "*")
title("Errore Energia meccanica, con N=400")
legend("eulero", "rk")
grid on
hold on
xlabel("t")
ylabel("E_meccanica")
