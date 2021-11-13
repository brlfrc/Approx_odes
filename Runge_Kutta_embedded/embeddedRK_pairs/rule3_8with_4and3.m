function [c,A,b,B,p,P]=rule3_8with_4and3

% embedded RK pair rule 3/8 and order 3 embedded

c=[0 1/3 2/3 1 1]';
A=[0 0 0 0 0;
   1/3 0 0 0 0;
   -1/3 1 0 0 0;
   1 -1 1 0 0;
   1/8 3/8 3/8 1/8 0];
b=[1/8 3/8 3/8 1/8 0]';
B=[1/12 1/2 1/4 0 1/6]';
p=4;
P=3;
