function [yy,nstep, nrej, nevals, Hr, Ha,STIMA,tt] = RKembedded(fun,t0,tf,y0,tableau,TOL)
[c,A,b,B,P,Q] = feval(tableau); %p,Q sono gli oridini dei due metodi
d = length(y0);
s = length(c); %numero di tappe
nB = length(B);
nb = length(b);
nstep = 0;
cost = 0;
nrej = 0;
h = max(TOL^(1/(P+1)),0.25);
tt = t0; 
t = t0;
K = zeros(d,s);
yy = y0;
yaux = y0;
yP = 0;
yQ = 0;
Ha = [1,1];
Hr = [1,1];
STIMA = [0,0];
count=0;
while(t<tf)
    if(t+h>tf)  %aggiusto il passo nel caso di uscita dal range
        h = tf-t;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %faccio i due Rk
    K(:,1) = fun(t,yaux);
    for i = 2:s
        somma = 0;
        for j = 1:i-1
            somma = somma + A(i,j)*K(:,j);
        end
        K(:,i) = fun(t+c(i)*h,yaux + h*somma);
    end
    cost = cost +s;
    sum = 0;
    for i=1:nb
        sum = sum + b(i)*K(:,i);
    end
    sum1 = 0;
    for i=1:nB
        sum1 = sum1 + B(i)*K(:,i);
    end
    yP = yaux + h*sum;
    yQ = yaux + h*sum1;
    
    est = norm(yP-yQ,'inf'); %sarebbe l'epsilon degli appunti
    if(est<TOL || count > 100) %accetto il passo
        Ha = [Ha;[t,h]];
        t = t+h;  %aumento di passo premipostato 
        nstep = nstep +1;
        yaux = yP; 
        yy = [yy,yP];
        tt = [tt,t];
        STIMA =[STIMA; [t, est]];
        count=0;
    else 
        nrej = nrej +1;
        Hr = [Hr;[t,h]];
        count = count +1;
    end

    hnew = h*0.75*max(0.4,(TOL/(est+eps))^(1/(min(P,Q)+1))); %vedi foglio per valore effettivo da inserire al posto di uno
    h = hnew;
end
Ha = Ha(2:end,:);
Hr = Hr(2:end,:);
nevals=cost;
end







% function [yy,nstep, nrej, nevals, H_r, H_a,STIMA,tt]= RKembedded (f,t0,tf,y0,tableau,TOL)
% 
% [c,A,b,B,p,P] = feval(tableau);
% H_a=[0,0]; % Passi accettati
% H_r=[0,0]; % Passi rifiutati
% STIMA =[0,0];
% tt=zeros(1, 1);
% tt(1)=t0;
% t_aux=t0;
% 
% n=length(y0);
% yy=zeros(n, 1);
% yy(:,1)=y0;
% y_aux=y0;
% 
% cost=0;
% m=0;
% nstep=0;
% nrej=0;
% h=max(TOL^(1/(p+1)),0.25);
% 
% while t_aux<tf
%     count=0;
%     
%     [y_aux0, y_aux1,cost_onestep,y_step]= RK_onestep (y_aux, t_aux,h, n, f, c,A,b,B);
%     est=max(abs(y_aux1-y_aux0));
%     cost= cost + cost_onestep;
% 
%     while (count<1000 && est>TOL)
%         h=h*0.75*max(0.4,min(1.5,(TOL/(est+0.01))^(1/(min(P,p)+1))));
%         
%         [y_aux0, y_aux1,cost_onestep,y_step]= RK_onestep (y_aux, t_aux,h, n, f, c,A,b,B);
% 
%         cost= cost + cost_onestep;
%         count=count+1;
%         nrej= nrej+1;
%         H_r = [H_r;[t_aux,h]];
%      
%         est=max(abs(y_aux1-y_aux0));
%     end
% 
%     if (t_aux+h>tf)
%         h=tf-t_aux;
%         [~, y_aux1,cost_onestep,y_step]= RK_onestep (y_aux, t_aux,h, n, f, c,A,b,B);
%         cost= cost + cost_onestep;
%     end
% 
%     y_aux= y_aux+h*y_step;
%     yy(:,m+1)=y_aux;
% 
%     t_aux=t_aux+h;
%     tt(m+1)=t_aux;
%     m=m+1;
% 
%     nstep=nstep+1;
%     H_a = [H_a;[t_aux,h]];
%     STIMA =[STIMA; [t_aux, max(abs(y_aux1-y_aux))]];
% end
% 
% nevals=cost;
% nrej=nrej-nstep;
% 
% end
% 
% function [y_aux0, y_aux1, cost,y_step]= RK_onestep(y_aux, t_aux, h, n, f, c,A,b,B)
% 
% s= size(A,1);
% K=zeros(n, s);
% cost=0;
% temp=0;
% 
% for i=1:s
%     for j=1:i
%         temp = temp + A(i,j)*K(:,j);
%     end
%     
%     K(:,i) = f(t_aux + c(i)*h, y_aux + temp*h);
%     cost=cost+1;
%     temp=0;
% end
% 
% temp1=0;
% for i=1:s
%         temp = temp +  b(i)*K(:,i);
%         temp1 = temp1 + B(i)*K(:,i);
% end
%     
% y_aux0 = y_aux+h*temp;
% y_aux1 = y_aux+h*temp1;
% 
% y_step=temp;
% end
