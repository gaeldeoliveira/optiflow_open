
function y = odeinit_blasius(x)
% initial guess for the solution
% Here we will use the classical polynomial approximation
%   f(eta) = 2*eta - eta^2
% Implying:
%   g(eta) = eta^2 - 1/3 * eta^3
% (divide by two because we also divided our equation!)

% Recast x into suitable orientation for vectorized evaluation
x = transpose(x(:));
x_in  = x(x<=1);
x_out = x(x>1);
% Compute g and its derivatives for initial solution!
% Case in which eta is below 1
y1_in = x_in.^2 - 1/3 * x_in.^3;
y2_in = 2*x_in  - x_in.^2;
y3_in = 2    - 2*x_in;

% For eta above 1
y1_out = 2/3 + (x_out-1);
y2_out = ones(size(x_out));
y3_out = zeros(size(x_out));

% Join all this into a single array!
y1 = [y1_in , y1_out];
y2 = [y2_in , y2_out];
y3 = [y3_in , y3_out];

% (divide by two because we also divided our equation!)
y = [y1 ; y2 ; y3];
end