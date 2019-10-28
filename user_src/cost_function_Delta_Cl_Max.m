function sc = cost_function_Delta_Cl_Max(~ , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors

%% Input extraction
% Extract parameters from parameters structure
    % No parameters for this case

% Extract aerodynamic polar from standard input form
result_soiled  = experiments_results{1};
ap_soiled = result_soiled.ap;

result_clean   = experiments_results{2};
ap_clean = result_clean.ap;

% Check that polars are valid
if and(strcmp(class(ap_clean) , 'aerodynamic_polar') , strcmp(class(ap_soiled) , 'aerodynamic_polar')) %#ok<STISA>
    % If polar is valid
    sc = ap_clean.cl_max_local_limited - ap_soiled.cl_max_local_limited;
else
    % If polar is invalid!
    sc = 0;
end