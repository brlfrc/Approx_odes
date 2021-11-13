function [c,A,b,B,p,P]=RRK_2and3

% embedded RK pair , order 2 and 3, from Iserles-book p. 84

c=[0 2/3 2/3]';
A=[0 0 0 ;
   2/3 0 0 ;
   0 2/3 0];
b=[1/4 3/4 0]';
B=[1/4 3/8 3/8]';
p=2;
P=3;
