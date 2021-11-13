function [c,A,b,B,p,P]=HIHA55_5and4

% embedded RK pair  Higham & Hall order 5 and 4 , from  Hairer-II, p. 27
% has FSAL property

c=[0 2/9 1/3 1/2 3/5 1 1]';
b=[1/12 0 27/32 -4/3 125/96 5/48 0]';
A=[zeros(1,7);
    2/9 zeros(1,6);
    1/12 1/4 zeros(1,5);
    1/8 0 3/8 zeros(1,4);
    91/500 -27/100 78/125 8/125 zeros(1,3);
    -11/20 27/20 12/5 -36/5 5 0 0;
    b'];
B=[2/15 0 27/80 -2/15 25/48 1/24 1/10]';
p=5;
P=4;
