function [c,A,b,B,p,P]=dopri5_5and4

% par  DoPri de orden 5 y 4 , tomado de  Hairer-I, p. 178
% o de Lambert, p. 186
% es FSAL

c=[0 1/5 3/10 4/5 8/9 1 1]';
b=[35/384 0 500/1113 125/192 -2187/6784 11/84 0]';
A=[zeros(1,7);
    1/5 zeros(1,6);
    3/40 9/40 zeros(1,5);
    44/45 -56/15 32/9 zeros(1,4);
    19372/6561 -25360/2187 64448/6561 -212/729 zeros(1,3);
    9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
    b'];
B=[5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]';
p=5;
P=4;
