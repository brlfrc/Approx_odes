clear all
close all

N_molecole=4;
q=zeros (2,N_molecole);
p=zeros (2,N_molecole);
q(:,1)=[1.45,0];
p(:,1)=[0,0];
q(:,2)=[1.3,1];
p(:,2)=[0,0];
q(:,3)=[-1,1.3];
p(:,3)=[0,0];
q(:,4)=[+1,-1];
p(:,4)=[0,0];
z0= convert(q,p);

t0=0;
tf=0.4;
h=[0.01,0.001,0.0001,0.00001];
x_box=3;
y_box=3;


U0=0.7*abs(U_potenziale_box(q,x_box,y_box));
fun=@(t,z) system_force_box (z,x_box,y_box,U0);

for j=1:length(h)
    
    [tt, yy,nevals_rk]= RKclassico (fun, t0, tf, h(j), z0, Verner6);
    yy=yy';
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

    diff(j)=abs(U(end)+K_media(end)-U(1));
end

figure
loglog(h,diff, "*")
xlabel("h")
ylabel("E(t_f)")
grid on

function z0=convert(q,p)
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