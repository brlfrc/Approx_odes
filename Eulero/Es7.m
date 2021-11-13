clear all
close all

% y'=[y2; -y1], y(0)=y0

i=1;

for N=[25,50,100,200,400]
    t0=0;
    tf1=10;
    tf2=40;
    h1=tf1/N;
    h2=tf2/N;
    
    y0=[pi/4; -pi/4];
    
    f=@(t,y) [y(2); -y(1)];

    [yy1,~,tt1]= euler_esplicito (f, t0, tf1, y0, h1);
    [yy2,~,tt2]= euler_esplicito (f, t0, tf2, y0, h2);
    [t_ode45,y_ode45] = ode45(f, [t0, tf2], y0);
    
    norm_10(i)=norm(yy1);
    i=i+1;

    figure(1)
    plot(tt1,yy1(:,1),"*")
    grid on
    hold on

    figure(2)
    plot(tt2,yy2(:,1),"*")
    grid on
    hold on
end
figure(1)
legend ("N= 25", "N= 50", "N= 100","N= 200", "N= 400")

figure(2)
legend ("N= 25", "N= 50", "N= 100","N= 200", "N= 400")

figure (3)
plot(yy2(:,1),yy2(:,2),"*", y_ode45(:,1), y_ode45(:,2), "*")
legend("Traiettoria", "traiettoria esatta")
grid on

figure (4)
title("Errore in norm_2 al variate di N (quindi h)")
loglog([25,50,100,200,400], norm_10,"*")
legend("errore")
xlabel("N")
ylabel("Errore")
grid on

%Non mi esce moto circolare