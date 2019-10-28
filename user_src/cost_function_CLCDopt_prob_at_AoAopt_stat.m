function sc = cost_function_CLCDopt_prob_at_AoAopt_stat(parameters , experiments_results)
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

try
    if strcmp(class(ap) , 'aerodynamic_polar') %#ok<STISA>
        % If polar is valid
        % Get range of angles of attack
        alpha_range = ap.alpha_range;
        % Get static L/D curve over it
        ld_alpha = ap.cl_alpha(alpha_range) ./ ap.cd_alpha(alpha_range);
        % Find maximum static L/D
        [ld_max_stat, i_ld_max_stat] = max(ld_alpha);
        % Find alpha opt of static L/D
        alpha_opt_stat = alpha_range(i_ld_max_stat);
        
        % Get expected L/D curve through convolution
        [alpha_bar_range, ld_expected] =   expected_polar_convolution(  sigma, n_sigma, n_alpha_til, alpha_range, ld_alpha);
        % Get expected L/D at max L/D point (from static case)
        ld_expected_at_alpha_opt_stat = interp1(alpha_bar_range, ld_expected, alpha_opt_stat);
        % Display
        disp(['Max L/D      = ' , num2str(ld_max_stat)])
        disp(['Expected L/D = ' , num2str(ld_expected_at_alpha_opt_stat)]);
        
        % Check that polar was converged enough
        if not(isnan(ld_expected_at_alpha_opt_stat))
            % Convolution was succesful
            sc = ld_expected_at_alpha_opt_stat;
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
catch
    % Something failed
    sc = 0;
    disp(['Something failed in: ', mfilename()]);
    disp(['x = '  , num2str(experiment_result.x)]);
end

% Make check to verify that sides do not cross
valid_shape = valid_shape_test(experiment_result);
if not(valid_shape)
   sc = 0; 
end


% Make additional check to verify that L/D did not enter into unreasonable
% branch