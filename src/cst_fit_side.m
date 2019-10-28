function [cst_parameters_vector , z_te , cst_parameters_vector_non_linear , cst_parameters_vector_full_non_linear] = cst_fit_side(cst_parametrization, tx_coordinates , tz_coordinates, a0)
%[ CST_PARAMETERS ] = CST_FIT_SIDE( TX , TZ, ORDER)
% Receives a set of adimensional airfoil side coordinates (tx and tz) and
% returns corresponding CST shape parameters
% 
% Inputs :
%   cst_parametrization - object of the parametrization class
%   including all necessary preprocessed data and methods for efficient
%   parametrization
%   tx  - a monotonic, crescent vector ordered in accordance with tz
%   tz  - is the adimensional thikness coordinate measured perpendicular to
%         the x axis (not the mean camber line of the airfoil)

%% Establish implicit interpolation functions
% This is what allows to deal with arbitrary order of vector. Could
% eventually lead to relaxation of monotonicity constraint but it is not a
% good practice to do so

% Currently working with default 'linear' interpolation. Not clear wether
% 'spline' or shape_preserving cubic 'pchip' would be better options
% Deprecate Original code to support loss of functionality in interp1
% pp = interp1(tx_coordinates, tz_coordinates, 'linear' , 'pp');
% Switch to pchip for R2013A onwards (LTS should become griddedinterpolant)
pp = pchip(tx_coordinates, tz_coordinates);

tz = @(tx) ppval(pp , tx);

% Interpolation function for derivative
pp_d1 = fnder(pp);
dz_dt = @(tx) ppval(pp_d1 , tx);


%% Initialize Parameter vectors

cst_parameters_vector = zeros(1, cst_parametrization.order);
cst_parameters_vector(1) = a0;


%% Obtain trailing edge data
% Determine Trailing Edge Bluntness
z_te = tz(1);

% Determine Boat tail angle
% Equals derivative dz_dt minus the trailing edge bluntness due to the
% contribution of the bluntness linear correction to the derivative of the
% final function
cst_parameters_vector(end) = - dz_dt(1); %- z_te;

% disp('z_te')
% disp(z_te)
% Clarika - Bien M'erit'e

%% Construct Point Vector
% With Correction for linear term added for blunt trailing edge
tz_fit_points_inner = tz(cst_parametrization.tx_positions_inner) - z_te * cst_parametrization.tx_positions_inner;

%% Construct RHS of linear system for solution
% disp(tz_fit_points_inner)
% disp(cst_parameters_vector(1) * cst_parametrization.sc_a0_line)
% disp(cst_parameters_vector(end) * cst_parametrization.sc_an_line)
RHS = transpose(tz_fit_points_inner) - cst_parameters_vector(1) * cst_parametrization.sc_a0_line - cst_parameters_vector(end) * cst_parametrization.sc_an_line;

%% Obtain Solution with precomputed inverse
cst_parameters_vector(2:(end-1)) = cst_parametrization.sc_submatrix_inverse * RHS;


%% Make Non-linear least squares fit
%  Now Make Non-linear least squares fit for inner points, based on
%  linear fit results

tz_coordinates_fit_nonlinear = tz_coordinates - z_te * tx_coordinates;

f = cst_parametrization.make_fit_object(cst_parameters_vector);

[c, cr] = fit(tx_coordinates, tz_coordinates_fit_nonlinear,f);


%% Now extract coefficients from fit structure f and fill
% cst_parameters_vector_non_linear vector

% Initialize vector
cst_parameters_vector_non_linear = zeros(size(cst_parameters_vector));
% Leading edge radius is the same
cst_parameters_vector_non_linear(1) = cst_parameters_vector(1);
% Trailing edge derivative is also the same
cst_parameters_vector_non_linear(end) = cst_parameters_vector(end);

% Fill vector (syntax is unclear but allows general degree fit)
for n_poly = 2:(length(cst_parameters_vector)-1)
    cst_parameters_vector_non_linear(n_poly) = eval(['c.a' , num2str(n_poly)]);
end


% Now Make fullly non-linear least squares fit

tz_coordinates_fit_nonlinear = tz_coordinates - z_te * tx_coordinates;

f = cst_parametrization.make_fit_object_full_non_linear(cst_parameters_vector);

[c, cr] = fit(tx_coordinates, tz_coordinates_fit_nonlinear,f);

cst_parameters_vector_full_non_linear = zeros(size(cst_parameters_vector));

for n_poly = 1:length(cst_parameters_vector)
    cst_parameters_vector_full_non_linear(n_poly) = eval(['c.a' , num2str(n_poly)]);
end

% Done!!!
% Les proverbes - Pierre Perret



end

