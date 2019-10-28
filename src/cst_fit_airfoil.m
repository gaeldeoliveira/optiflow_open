function [parameters_upper , z_te_upper , parameters_lower , z_te_lower , parameters_upper_non_linear , parameters_lower_non_linear , parameters_upper_vector_full_non_linear , parameters_lower_vector_full_non_linear] = cst_fit_airfoil(cst_parametrization_upper , cst_parametrization_lower , tx_coordinates , tz_coordinates)
%CST_FIT(TX_COORDINATES, TZ_COORDINATES , [upper_order] , [lower_order])
%This function receives a set of airfoil coordinates and returns a set of
%CST parameters.
%   Compulsory Arguments: 
%       TX_COORDINATES - Normalized chordwise coordinates in the order
%         1->0->1, that is, from trailing edge to leading edge by the upper
%         side and comming back to trailing edge by the lower side
%       TZ_COORDINATES - Normalized Z coordinates in cartesian space (not
%         from mean camber line)
%   Optional arguments:
%       [upper_order] - Order of Bersntein Polynomials from upper side
%       [lower_order] - Order of Bersntein Polynomials from lower side
%
%   When no optional arguments are supplied the order is set to 5 on both
%   sides by default. If [upper_order] is specified then [lower_order]
%   should be specified as well.


%% Separate Upper and Lower Sides
% Find leading edge index
[~ , min_index] = min(tx_coordinates);

% Generate Upper and Lower tx coordinate vectors
tx_lower_side = tx_coordinates(min_index:end);
tx_upper_side = flipud(tx_coordinates(1:min_index)); % Reorder upper side to avoid confusions further down

% Same for tz coordinate vectors
tz_lower_side = - tz_coordinates(min_index:end);
tz_upper_side = flipud(tz_coordinates(1:min_index)); % Reorder upper side to avoid confusions further down

% Now extract leading edge data for radius of curvature determination
le_points_each_side = 12;
tx_le =  tx_coordinates((min_index-le_points_each_side) : (min_index+le_points_each_side));
tz_le =  tz_coordinates((min_index-le_points_each_side) : (min_index+le_points_each_side));



%% Extract Leading Edge Radius
% Make nonlinear least squares fit to obtain close approximation to le
% radius, with 12 points on each side, as in NASA TM 4906
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',0,...
               'Upper',20, ...
               'Startpoint',1);
f = fittype('b*z^2' , 'independent', {'z'} , 'options' , s);

[c,~] = fit(tz_le, tx_le,f);

a0_least_squares = 1/sqrt(c.b);

%% Fit Each Side

[parameters_upper , z_te_upper , parameters_upper_non_linear , parameters_upper_vector_full_non_linear] = cst_fit_side(cst_parametrization_upper , tx_upper_side , tz_upper_side , a0_least_squares);
[parameters_lower , z_te_lower , parameters_lower_non_linear , parameters_lower_vector_full_non_linear] = cst_fit_side(cst_parametrization_lower , tx_lower_side , tz_lower_side , a0_least_squares);


% Done !!!
% Blanche - Pierre Perret
end

