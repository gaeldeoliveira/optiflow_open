function [ fun_integral ] = riemann_integral(fun, y)
%SIMPLE_INTEGRAL performs a simple riemann integral numerically
%        int(fun, dy, y(1),y(end))
% Where:
%   fun     is a function handle to integrand (vectorized)
%   y       is a vector of intervals over which to conduct the integration

% Determine Spacing of intervals
y_spacing = y(2:end) -y(1:(end)-1);
% Determine Centers of intervals
y_centers = y(1:(end)-1) + y_spacing/2;
% Compute Contribution of each interval (evaluated in the middle)
fun_integral_steps = fun(y_centers) .* y_spacing;
% Sum up to obtain integral value
fun_integral = sum(fun_integral_steps);

end

