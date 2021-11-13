function [AA]=Heun3

% Tableu Heun  order 3 as in Hairer-I, p. 135

c=[0 1/3 2/3]';
b=[1/4 0 3/4]';
A=[0 0 0 ;
   1/3 0 0 ;
   0 2/3 0];

n=size(A,1)+1;
AA=zeros(n,n);
AA(1:n-1,1)=c;
AA(n,2:n)=b';
AA(1:n-1,2:n)=A;