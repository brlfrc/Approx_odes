function [Z,nevals] = newton(F, dF, Z0, atol, maxiter,Yn,h,t_aux)

    Z = Z0;
    theta = 0.5;
    i = 1;
    nevals=0;
    Delta = - dF(Z,h,t_aux) \ F(Z,Yn,h,t_aux);
    Z = Z + Delta;

    while norm(Delta)*theta/(1-theta) > atol && i < maxiter
    
        Delta_old = Delta;
        Delta = - dF(Z,h,t_aux) \ F(Z,Yn,h,t_aux);
        Z = Z - Delta;
        theta = norm(Delta)/norm(Delta_old);
        i = i+1;
        nevals=nevals+2;
    end

end