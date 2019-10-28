function sc = cost_function_CLCDopt_prob(parameters , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors

%% Input extraction
% Extract parameters from parameters structure
n_alpha_til = parameters.n_alpha_til;       % n_alpha_til = 64;          % [int. ] - number of elements convolution integral kernel
n_sigma     = parameters.n_sigma    ;       % n_sigma     = 1.64485;     % [adim.] - Width of integration range (in sigmas, affects pdf/cdf coverage)

% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{1};
ap = experiment_result.ap;

%% Extract last coordinate (equal to turbulence intensity)
x = experiment_result.x;
sigma       = x(end);      % [deg  ] - Standard deviation of angle of attack perturbations


%% Now run usual clcd at cl cost function


if strcmp(class(ap) , 'aerodynamic_polar')
    % If polar is valid
    % Get range of angles of attack
    alpha_range = ap.alpha_range;
    % Get static L/D curve over it
    ld_alpha = ap.cl_alpha(alpha_range) ./ ap.cd_alpha(alpha_range);
    % Get expected L/D curve through convolution
    [~, ld_expected] =   expected_polar_convolution(  sigma, n_sigma, n_alpha_til, alpha_range, ld_alpha);
    % Get expected L/D curve through convolution
    ld_expected_max = max(ld_expected);
    
    % Check that polar was converged enough 
    if not(isnan(ld_expected_max))
        % Convolution was succesful
        sc = ld_expected_max;
    else
        % Convolution failed, translate nan into 0
        sc = 0;
    end
   
    % Rule out "ultra-high performance" unphysical solutions (L/D > 300)
    if max(abs(ap.raw_data.cl ./ ap.raw_data.cd)) > 300
        sc = 0;
    end
else
    % If polar is invalid!
    sc = 0;
end

% Make check to verify that sides do not cross
valid_shape = valid_shape_test(experiment_result);
if not(valid_shape)
   sc = 0; 
end

% Make additional check to verify that L/D did not enter into unreasonable
% branch