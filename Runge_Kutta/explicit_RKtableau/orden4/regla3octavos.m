function AA =regla3octavos

% Tableu  3/8 rule order 4 as in Hairer-I, p. 138

c=[0 1/3 2/3 1]';
b=[1/8 3/8 3/8 1/8]';
A=[0 0 0 0;
   1/3 0 0 0;
   -1/3 1 0 0;
   1 -1 1 0];

n=size(A,1)+1;
AA=zeros(n,n);
AA(1:n-1,1)=c;
AA(n,2:n)=b';
AA(1:n-1,2:n)=A;