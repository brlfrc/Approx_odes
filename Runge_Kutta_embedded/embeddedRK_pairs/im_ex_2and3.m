function [c,A,b,B,p,P]=im_ex_2and3

% embedded RK pair RK, order 2 and 3, from Iserles book p. 89, ej. 5.5
% the order 2 is implicit
% the order 3 is explicit

c=[0 1 1/2]';
A=[0 0 0 ;
   1/2 1/2 0 ;
   3/2 -1 0];
b=[1/2 1/2 0]';
B=[1/6 2/3 1/6]';
p=2;
P=3;
