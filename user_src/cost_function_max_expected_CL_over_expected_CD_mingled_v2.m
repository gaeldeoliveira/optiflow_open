function sc = cost_function_max_expected_CL_over_expected_CD_mingled_v2(parameters , experiments_results)
% The cost_function_max_expected_CL_over_expected_CD_mingle objective function
% mingles results from cost_function_max_expected_CL_over_expected_CD
% obtained from experiment number 1 and 2, usually corresponding to free
% and forced transition conditions

% Define Mingling factor (0 = fully free, 1=fully tripped)
% mingle_factor = parameters.mingle_factor;

% Compute free transition objective value from polar
parameters_free = parameters;
parameters_free.n_experiment = 1;
[sc_free , alpha_design]= cost_function_max_expected_CL_over_expected_CD_v2(parameters_free , experiments_results);

% Compute, if results make any sense, compute roughness sensitivity
if not(isnan(alpha_design))
    % Extract aerodynamic polars
    % Free
    experiment_result_free   = experiments_results{1};
    ap_free                  = experiment_result_free.ap; 
    % Tripped
    experiment_result_trip   = experiments_results{2};
    ap_trip                  = experiment_result_trip.ap; 
    % Check that both ap objects are polars
    if and(strcmp(class(ap_free) , 'aerodynamic_polar'), strcmp(class(ap_trip) , 'aerodynamic_polar')) %#ok<STISA>
        % In that case, compute cl at alpha_design
        cl_design_free = ap_free.cl_alpha(alpha_design);
        cl_design_trip = ap_trip.cl_alpha(alpha_design);
        % Now check that cl is valid in both cases
        if and(isfinite(cl_design_free), isfinite(cl_design_trip))
            % In that case, compute cl difference
            delta_cl_trip = cl_design_free - cl_design_trip;
        else 
            delta_cl_trip = nan;
        end
    else
        delta_cl_trip = nan;
    end
else
    delta_cl_trip = nan;
end


% Now, check if delta_cl_trip is acceptable
if isnan(delta_cl_trip)
    % If not, return invalid flag value 
    sc = 0;
else
    % Otherwise
    if abs(delta_cl_trip) < 0.05
        % Return free performance index if roughness sensitivity meets
        % expectations
        sc = sc_free;
    elseif abs(delta_cl_trip) < 0.10
        % Roughness sensitivity is not yet as desired but still reasonable
        % Penalize cost function in substantial way
        sc = sc_free - 2000*abs(delta_cl_trip);
        % Remark that
        % sc_free is 100 to 300 range
        %   2000 * 0.05 = 100
        %   2000 * 0.10 = 200
    else
        % Roughness sensitivity is way off target, eliminate design
        sc = 0;
    end
end

% Make sure we never allow negative expected L/Ds due to penalization or
% other weird stuff
if sc < 0
    sc = 0;
end

end