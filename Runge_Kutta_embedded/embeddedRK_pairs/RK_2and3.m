function [c,A,b,B,p,P]=RK_2and3

% embedded RK pair, oders  2 and 3, from Iserles-book p. 89, ej. 5.4

c=[0 1/2 1]';
A=[0 0 0 ;
   1/2 0 0 ;
   -1 2 0];
b=[0 1 0]';
B=[1/6 2/3 1/6]';
p=2;
P=3;
