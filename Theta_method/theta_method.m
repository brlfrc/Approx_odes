function [ t, y ] = theta_method ( f, tspan, y0, n, theta )

%*****************************************************************************80
%
%% theta_method uses the theta method method + fsolve() to solve an ODE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 March 2021
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    function handle f: evaluates the right hand side of the ODE.  
%
%    real tspan(2): the starting and ending times.
%
%    real y0(m): the initial conditions. 
%
%    integer n: the number of steps.
%
%    real theta: the value of theta.
%    0 <= theta <= 1.
%
%  Output:
%
%    real t(n+1,1), y(n+1,m): the solution estimates.
%
  options = optimoptions ( 'fsolve', 'Display', 'off' );

  m = length ( y0 );
  t = zeros ( n + 1, 1 );
  y = zeros ( n + 1, m );

  dt = ( tspan(2) - tspan(1) ) / n;

  t(1,1) = tspan(1);
  y(1,:) = y0(:);

  for i = 1 : n

    to = t(i,1);
    yo = y(i,:);

    tn = to + dt; 
    yn = yo + dt * transpose ( f ( to, yo ) );

    yn = fsolve ( @(yn)theta_residual(f,to,yo,tn,yn,theta), yn, options );

    t(i+1,1) = tn;
    y(i+1,:) = yn;

  end

  return
end
function value = theta_residual ( f, to, yo, tn, yn, theta )

%*****************************************************************************80
%
%% theta_residual evaluates the theta method residual.
%
%  Discussion:
%
%    We are seeking a value YN defined by the implicit equation:
%
%      YN = YO + ( tn - to ) 
%        * ( theta * F ( TO, YO ) + ( 1.0 - theta ) * F ( TN, YN ) )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 March 2021
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    function handle f: evaluates the right hand side of the ODE.  
%
%    real to, yo: the old time and solution value.
%
%    real tn, yn: the new time and solution value.
%
%    real theta: the value of theta.
%    0 <= theta <= 1.
%
%  Output:
%
%    real value: the theta method residual.
%
  value = yn - yo - ( tn - to ) ...
    * transpose ( theta * f ( to, yo ) + ( 1.0 - theta ) * f ( tn, yn ) );

  return
end

