function AA =Runge3

% Tableu Runge3  order 3 as in  Hairer-I, p. 135

c=[0 1/2 1 1]';
b=[1/6 2/3 0 1/6]';
A=[0 0 0 0;
   1/2 0 0 0;
   0 1 0 0;
   0 0 1 0];

n=size(A,1)+1;
AA=zeros(n,n);
AA(1:n-1,1)=c;
AA(n,2:n)=b';
AA(1:n-1,2:n)=A;