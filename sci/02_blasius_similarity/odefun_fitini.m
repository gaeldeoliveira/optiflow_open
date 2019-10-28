
function y_prime = odefun_fitini(eta, y, lambda, t_p)
% Derivative function. (inspired on the mat4bvp example script!)
%   Fitini, la gentille équation!
% The sweet equation from the report!

% Compute omega with plasma Y weighting function! (imported for fitini <3,
% la gentille équation, qui donne également son nom a mon levain!)
omega  = wy_function_fitini(eta, t_p) ./ t_p;

% Now compute Fitini =)
y1_prime = y(2);
y2_prime = y(3);
y3_prime = - 0.5 * (y(1) * y(3) + 2 * omega * lambda * (1-lambda)) / (1-lambda).^2;

y_prime = [ y1_prime ; ...
            y2_prime ; ...
            y3_prime];

% x remains unused
end