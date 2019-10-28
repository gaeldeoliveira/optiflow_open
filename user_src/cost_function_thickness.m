function thickness = cost_function_thickness( ~ , experiments_results)
% Calculated weighted Cl/Cd of a polar at AOA/Weight combinations from
% alpha_i and w_i vectors

% Input extraction
% parameters can be ignored

% Extract aerodynamic polar from standard experiment input form
experiment_result   = experiments_results{1};
% parameters = experiment_result.x;
coordinates = experiment_result.coordinates;

tz = coordinates.tz;

%% Compute Thickness
% Projected Thickness
thickness = max(tz)-min(tz);

end