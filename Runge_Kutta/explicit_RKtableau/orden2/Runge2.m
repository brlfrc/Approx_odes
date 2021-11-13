function AA =Runge2

%  Runge order 2  Hairer-I, p. 135
% Called RK22 in Butcher, p. 87
c=[0 1/2]';
b=[0 1]';
A=[0 0;
   1/2 0];

n=size(A,1)+1;
AA=zeros(n,n);
AA(1:n-1,1)=c;
AA(n,2:n)=b';
AA(1:n-1,2:n)=A;
end