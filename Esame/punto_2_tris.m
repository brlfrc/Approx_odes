clear all
close all

t0=1;
tf=2;
a=0;
h=0.001;
fun=@(t,z) system_force(z);
fun1=@(t,z) system_force_modified(z);
for j=9:2:19
    a=a+1;
    N_molecole=(j+1)*(j+1);
    for i=0:j
        q(:,(i+1)+j*i:(i+1)+j*i+j)=[i+zeros(1,j+1);
            0:j];
    end
    
    p=zeros(2,N_molecole);
    z0= convert(q,p);

    tic
    [yy,nevals,tt, count1] = Storme_Verlet(fun, t0, tf, z0, h);
    time_normal(a)=count1;
    tic
    [yy,nevals,tt, count2] = Storme_Verlet(fun1, t0, tf, z0, h);
    time_modified(a)=count2;
end

figure
semilogy((10:2:20).^2,time_normal, "*",(10:2:20).^2,time_modified, "*")
xlabel("Numero di particelle")
ylabel("Valutazione forza")
grid on

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