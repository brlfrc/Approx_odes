clear all
close all

N_molecole=16;
x_box=4;
y_box=4;
q=position(x_box, y_box);
p=-0.5+rand(2,N_molecole);

figure()
hold on
for i=1:N_molecole
    quiver(q(1,i),q(2,i),p(1,i),p(2,i));
    plot(q(1,i),q(2,i),"*");
end
x = [2, -2, -2, 2, 2]; y = [2, 2, -2, -2, 2];
plot(x, y, 'b-', 'LineWidth', 3);


z0= convert(q,p);
t0=0;
tf=5;
h=0.0001;

K0=0;
for i=1:N_molecole
    K0=K0+0.5*(p(1,i).^2+p(2,i).^2);
end
U0=0.7*abs(U_potenziale_box(q,x_box,y_box)+K0);
fun=@(t,z) system_force_box (z,x_box,y_box,U0);
[tt, yy,nevals_rk]= RKclassico_box (fun, t0, tf, h, z0, Verner6);
yy=yy';

[K_media,U]=plot_Tra_EN (tt,yy, N_molecole,x_box,y_box,U0);

figure
semilogy(tt(1:end-1),abs(tt(1:end-1)-tt(2:end)),"*")
grid on
xlabel("t")
ylabel("h")

function [q]=position(x_box, y_box)
    a=1;
    for i= -x_box/2+0.1:1+0.1*rand():x_box/2-0.1
        for j= -y_box/2+0.1:1+0.1*rand():y_box/2-0.1
            q(:,a)=[i+0.1*rand(); j+0.1*rand()];
            a=a+1;
        end
    end
end

function z0= convert(q,p)
    i=1;
    for a=1:length(q)
        z0(i:i+1,1)=q(:,a);
        i=i+2;
    end
    for a=1:length(q)
        z0(i:i+1,1)=p(:,a);
        i=i+2;
    end
end

function [K_media,U]= plot_Tra_EN(tt,yy, N_molecole,x_box,y_box,U0)
    figure()
    x = [2, -2, -2, 2, 2]; y = [2, 2, -2, -2, 2];
    plot(x, y, 'b-', 'LineWidth', 3);
    for i=1:2:N_molecole*2
        hold on
        plot(yy(:,i),yy(:,i+1),"*");
    end
    grid on;
    
    for a=1:length(tt)
        K_media(a)=0;
        for i=N_molecole*2+1:2:N_molecole*4
           K_media(a)=K_media(a)+0.5*(yy(a,i).^2+yy(a,i+1).^2);
        end
    end

    for a=1:length(tt)
        b=0;
        for i=1:2:N_molecole*2
            b=b+1;
            q_tt(1,b)=yy(a,i);
            q_tt(2,b)=yy(a,i+1);
        end
        U(a)=U_potenziale_box(q_tt,x_box,y_box,U0);
    end
    figure()
    plot(tt,K_media,"*")
    hold on; grid on;
    plot(tt,U,"*")
    xlabel("t"); ylabel("E"); legend("Cinetica", "Potenziale");
    
    figure
    semilogy(tt,abs((K_media+U)-U(1)-K_media(1)), "*")
    grid on; xlabel("t"); ylabel("\DeltaE(t)");

    figure
    E=abs(U(1)+K_media(1));
    semilogy(tt,abs((K_media+U)-U(1)-K_media(1))/E, "*")
    grid on; xlabel("t"); ylabel("\DeltaE_%(t)");
end
