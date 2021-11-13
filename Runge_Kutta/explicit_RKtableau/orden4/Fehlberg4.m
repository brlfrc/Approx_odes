function AA =Fehlberg4

% Tableau Fehlberg  order 4 as in Hairer-I, p. 177

c=[0 1/4 3/8 12/13 1 1/2]';
A=[zeros(1,6);
    1/4 zeros(1,5);
    3/32 9/32 zeros(1,4);
    1932/2197 -7200/2197 7296/2197  zeros(1,3);
    439/216 -8 3680/513 -845/4104 zeros(1,2);
    -8/27 2 -3544/2565 1859/4104 -11/40 0];
b=[25/216 0 1408/2565 2197/4104 -1/5 0]';

n=size(A,1)+1;
AA=zeros(n,n);
AA(1:n-1,1)=c;
AA(n,2:n)=b';
AA(1:n-1,2:n)=A;