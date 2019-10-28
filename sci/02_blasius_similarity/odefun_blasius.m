
function dydx = odefun_blasius(x, y)
% Derivative function. (inspired on the mat4bvp example script!)
% The classical blasius equation reads
%       g'''+ 0.5*g*g''=0
% As described in the report, we recast it into a system of 3 non-linear
% first order ODEs

dydx = [ y(2) ; ...
    y(3) ; ...
    - 0.5 * y(1) * y(3)];

% Lambda and x remain unused
end