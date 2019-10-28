function sc = cost_function_max_expected_CL_over_expected_CD_mingled(parameters , experiments_results)
% The cost_function_max_expected_CL_over_expected_CD_mingle objective function
% mingles results from cost_function_max_expected_CL_over_expected_CD
% obtained from experiment number 1 and 2, usually corresponding to free
% and forced transition conditions

% Define Mingling factor (0 = fully free, 1=fully tripped)
mingle_factor = parameters.mingle_factor;

% Compute free transition objective value from polar
parameters_free = parameters;
parameters_free.n_experiment = 1;
sc_free = cost_function_max_expected_CL_over_expected_CD(parameters_free , experiments_results);
% Compute tripped transition objective value from polar
parameters_trip = parameters;
parameters_trip.n_experiment = 2;
sc_trip = cost_function_max_expected_CL_over_expected_CD(parameters_trip , experiments_results);
% Mingle cost function
sc      = (1 - mingle_factor) * sc_free + mingle_factor * sc_trip;

end