function [Z] = newton(F, dF, Z0, atol, maxiter,yn,h,tn,c,fun,Jac)

    Z = Z0;
    theta = 0.5;
    i = 1;
    Delta = - dF(Z,h,Jac) \ F(Z,h,tn,yn,c,fun);
    Z = Z + Delta;

    while norm(Delta)*theta/(1-theta) > atol && i < maxiter
    
        Delta_old = Delta;
        Delta = - dF(Z,h,Jac) \ F(Z,h,tn,yn,c,fun);
        Z = Z - Delta;
        theta = norm(Delta)/norm(Delta_old);
        i = i+1;
    end

end