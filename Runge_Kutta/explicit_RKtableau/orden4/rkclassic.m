function AA =rkclassic

% Tableau classic RK clasico de orden 4 como en Hairer-I, p. 138

c=[0 1/2 1/2 1]';
b=[1/6 1/3 1/3 1/6]';
A=[0 0 0 0;
   1/2 0 0 0;
   0 1/2 0 0;
   0 0 1 0];

n=size(A,1)+1;
AA=zeros(n,n);
AA(1:n-1,1)=c;
AA(n,2:n)=b';
AA(1:n-1,2:n)=A;