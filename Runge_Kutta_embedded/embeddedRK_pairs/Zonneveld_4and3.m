function [c,A,b,B,p,P]=Zonneveld_4and3

% embedded RK pair Zonneveld 4(3)
% as in Hairer-I, p. 167

c=[0 1/2 1/2 1 3/4]';
A=[0 0 0 0 0;
   1/2 0 0 0 0;
   0 1/2 0 0 0;
   0 0 1 0 0 ;
   5/32 7/32 13/32 -1/32 0];
b=[1/6 1/3 1/3 1/6 0]';
B=[-1/2 7/3 7/3 13/6 -16/3]';
p=4;
P=3;
