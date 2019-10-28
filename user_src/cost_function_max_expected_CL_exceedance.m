function sc = cost_function_max_expected_CL_exceedance(parameters , experiments_results)
%function sc = cost_function_max_expected_CL_over_expected_CD
% Used in probabilistic design (chap 7 of PhD thesis)
% Return relative stat cl difference between cl at alpha_design  and cl at alpha_design + n_sigma 
% alpha_design is the one that maximizes l_expected / d_expected

%% Preparameter extraction 
% (to select experiment_result if applicable)
if isfield(parameters, 'n_experiment')
    n_experiment = parameters.n_experiment;
else
    n_experiment = 1;
end

%% Experiment extraction
% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{n_experiment};
ap = experiment_result.ap; 

% Get sigma_std from parameters
std_alpha_deg = parameters.sigma;
n_sigma_gust  = parameters.n_sigma_gust;

% Define invalidity flag
sc_invalid_value = inf;

%% Process polar

% Check that this is really an aerodynamics polar
if strcmp(class(ap) , 'aerodynamic_polar') %#ok<STISA>
    % If it is, get range of converged angles of attack
    alpha_range       = ap.alpha_range();
    % Check its width
    alpha_range_width = max(alpha_range) - min(alpha_range);
    
    % Airfoil is considered invalid if converged range is smaller than 4 sigma
    if alpha_range_width < 6 * std_alpha_deg
        sc = sc_invalid_value;
    else
        % Make sanitized range of angles of attack
        alpha_range_sanitized = alpha_range(and(alpha_range>(min(alpha_range)+2.5*std_alpha_deg), alpha_range<(max(alpha_range)-2.5*std_alpha_deg)));
        
        % Get expected cl
        cl_expected = ap.field_alpha_prob('cl', alpha_range_sanitized, std_alpha_deg);
        % Get expected cd
        cd_expected = ap.field_alpha_prob('cd', alpha_range_sanitized, std_alpha_deg);
        % Compute ratio
        cl_expected_over_cd_expected = (cl_expected ./ cd_expected);
        % Compute maximum value of l_expected / cd_expected
        [~, i_alpha_range_sanitized] = max(cl_expected_over_cd_expected);
        % Determine corresponding design angle of attack
        alpha_design     = alpha_range_sanitized(i_alpha_range_sanitized);
        % Determine gust angle of attack
        alpha_gust       = alpha_design + n_sigma_gust * std_alpha_deg;
        % Determine stat cl at design angle of attack
        cl_stat_design   = ap.cl_alpha(alpha_design);
        % Determine max stat cl between alpha_design and alpha_gust
        alpha_gust_range = linspace(alpha_design, alpha_gust, 32);
        % Determine stat cl within gust range
        cl_gust_range     = ap.cl_alpha(alpha_gust_range);
        % Determine max stat cl within gust range
        cl_gust_max       = max(cl_gust_range);
        % Determine max cl design exceedance (0 is desired value)
        sc = cl_gust_max - cl_stat_design;
    end
    
    % Rule out "ultra-high performance" unphysical solutions (L/D > 300 in free transition)
    if max(abs(ap.raw_data.cl ./ ap.raw_data.cd)) > 300
        sc = sc_invalid_value;
    end
    % Rule out "ultra-high performance" unphysical solutions (L/D > 200 in forced transition)
    if and(max(abs(ap.raw_data.cl ./ ap.raw_data.cd)) > 160, n_experiment == 2)
        sc = sc_invalid_value;
    end
    
else
    sc = sc_invalid_value;
end

% Make check to verify that trailing edge is thick enough
[pn75, pn90] = low_thickness_exceedance(parameters, experiment_result);
if or(pn75 > 0, pn90 > 0 )
    % Limit pn factors to 0 (shape crossing handled separately)
    pn75 = min(pn75, 1);   % 0 means no exceedance, 1 means full exceedance
    pn90 = min(pn90, 1);   % 0 means no exceedance, 1 means full exceedance
    % Damp out expected L over D
    % sc = sc * (1 - pn75).^2 * (1 - pn90).^2;
    % Artificially augment expected Cl_exceedance
    sc = sc + (1 / ((1 - pn75).^2 * (1 - pn90).^2) - 1);
end

% Eliminate overly favourable cases
if sc < 0.02
    sc = sc_invalid_value;
end

%% Make check to verify that sides do not cross
valid_shape = valid_shape_test(experiment_result);
if not(valid_shape)
   sc = sc_invalid_value;
end


% Make additional check to verify that L/D did not enter into unreasonable
% branch