function [c,A,b,B,p,P]=EulerHeun

% embedded RK pair with Euler and ''improved- Euler''
c=[0 1]';
A=[0 0;
   1 0];
b=[1 0]';
B=[1/2 1/2]';
p=1;
P=2;
