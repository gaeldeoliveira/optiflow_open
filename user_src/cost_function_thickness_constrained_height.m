function thickness = cost_function_thickness_constrained_height( parameters , experiments_results)
% Calculated weighted Cl/Cd of a polar at AOA/Weight combinations from
% alpha_i and w_i vectors

%% Extract Inputs
% Extract parameters from parameters structure
tx_no_check   = parameters.tx_no_check  ;
min_thickness = parameters.min_thickness;
max_thickness = parameters.max_thickness;

% Extract experiment result from cell array (even if it is a one entry
% array!)
experiment_result   = experiments_results{1};
% Get the airfoil coordinates out of experiment result
coordinates = experiment_result.coordinates;

tz = coordinates.tz;

%% Compute Thickness
% Projected Thickness
thickness = max(tz)-min(tz);

%% Compute Penalty Factor for Minimum Building Height
penalty_factor = penalty_factor_building_height_violation(tx_no_check, min_thickness, coordinates);

%% Compute Penalty Cost for Exceeding Maximum Building Height
if thickness > max_thickness
    penalty_cost_max_thickness = thickness - max_thickness;
else
    penalty_cost_max_thickness = 0;
end

%% Apply penalty factor and cost to thickness 
%(make airfoil look thinner than it is if a constraint is violated)
thickness      = penalty_factor * thickness - 2 * penalty_cost_max_thickness;

end