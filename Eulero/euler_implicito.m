function [t,y] = euler_back(t0,y0,t_end,h,fcn,tol)
%----------------------------------------------------
% MATLAB CODE TAKEN FROM
% K. Atkinson, W. Han e D. Steward,
% Numerical Solution of Ordinary Differential Equations,
% programs for the book, John Wiley and Sons, 2009,
% http://www.math.uiowa.edu/NumericalAnalysisODE/
%----------------------------------------------------
% function [t,y] = euler_back(t0,y0,t_end,h,fcn,tol)
%----------------------------------------------------
% Solve the initial value problem
% y’ = f(t,y), t0 <= t <= b, y(t0)=y0
% Use the backward Euler method with a stepsize
% of h.
% The user must supply an m-file to define the
% derivative f, with some name,
% say ’deriv.m’, and a first line of the form
% function ans=deriv(t,y)
% tol is the user supplied bound on the difference
% between successive values of the backward Euler
% iteration.
% A sample call would be
% [t,z]=euler_back(t0,z0,b,delta,’deriv’,1e-6)
%----------------------------------------------------
% Output:
% The routine euler_back will return two vectors,
% t and y.
% The vector t will contain the node points
% t(1)=t0, t(j)=t0+(j-1)*h, j=1,2,...,N
% with
% t(N) <= t_end, t(N)+h > t_end
% The vector y will contain the estimates of the
% solution Y at the node points in t.
%----------------------------------------------------
    % Initialize.
    n = fix((t_end-t0)/h) + 1;
    t = linspace(t0,t0+(n-1)*h,n)';
    y = zeros(n,1);
    y(1) = y0;
    % advancing
    for i=2:n
        % forward Euler estimate
        yt1 = y(i-1) + h*feval(fcn,t(i-1),y(i-1));
        % one-point iteration
        count = 0;
        diff = 1;
        while diff > tol && count < 10
            yt2 = y(i-1) + h*feval(fcn,t(i),yt1);
            diff = abs(yt2-yt1);
            yt1 = yt2;
            count = count +1;
        end
        if count >= 10
            disp("Not converging after 10 steps at t = "+num2str(t(i)))
        end
        y(i) = yt2;
    end
end % euler_back