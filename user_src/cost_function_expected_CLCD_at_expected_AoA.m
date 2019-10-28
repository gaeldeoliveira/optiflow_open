function sc = cost_function_expected_CLCD_at_expected_AoA(parameters , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors


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

%% Parameter extraction
% Extract parameters from parameters structure
alpha_bar   = parameters.alpha_bar  ;       % alpha_bar   = 7 ;          % [deg. ] - prescribed average angle of attack
n_alpha_til = parameters.n_alpha_til;       % n_alpha_til = 64;          % [int. ] - number of elements convolution integral kernel
n_sigma     = parameters.n_sigma    ;       % n_sigma     = 1.64485;     % [adim.] - Width of integration range (in sigmas, affects pdf/cdf coverage)

% Get sigma                                 % sigma       = 2 ;          % [deg  ] - Standard deviation of angle of attack perturbations
if isfield(parameters, 'sigma')
    % From parameters structure if applicable
    sigma   = parameters.sigma;
else
    % From genotype if not specified in parameters structure
    % (corresponds to last dummy parameter by convention)
    sigma   = experiment_result.x(end);
end


%% Now run usual cost function computations
try
    if strcmp(class(ap) , 'aerodynamic_polar') %#ok<STISA>
        % If polar is valid
        % Get range of angles of attack
        alpha_range = ap.alpha_range;
        % Get static L/D curve over it
        ld_alpha = ap.cl_alpha(alpha_range) ./ ap.cd_alpha(alpha_range);
        % Find maximum static L/D
        [ld_max_stat, i_ld_max_stat] = max(ld_alpha);

        % Get expected L/D curve through convolution
        [alpha_bar_range, ld_expected] =   expected_polar_convolution(  sigma, n_sigma, n_alpha_til, alpha_range, ld_alpha);
        % Get expected L/D at max L/D point (from static case)
        ld_expected_at_alpha_bar = interp1(alpha_bar_range, ld_expected, alpha_bar);
        % Display
        disp(['Max L/D expt = ' , num2str(max(ld_expected)         , 4 ), '       with sigma = ', num2str(sigma)                        ])
        disp(['Max L/D stat = ' , num2str(ld_max_stat              , 4) , '     at alpha_opt = ', num2str(alpha_range(i_ld_max_stat), 3)]);
        disp(['Expected L/D = ' , num2str(ld_expected_at_alpha_bar , 4) , '     at alpha_bar = ', num2str(alpha_bar                 , 3)]);
        
        % Check that polar was converged enough
        if not(isnan(ld_expected_at_alpha_bar))
            % Convolution was succesful
            sc = ld_expected_at_alpha_bar;
        else
            % Convolution failed, translate nan into 0
            sc = 0;
        end
        
        % Rule out "ultra-high performance" unphysical solutions (L/D > 350)
        if max(abs(ap.raw_data.cl ./ ap.raw_data.cd)) > 350
            sc = 0;
        end
    else
        % If polar is invalid!
        sc = 0;
    end
catch
    % Something failed
    sc = 0;
    lg.err(ap, ['Something failed in: ', mfilename()]);
    lg.err(ap, ['x = '  , num2str(experiment_result.x)]);
end

% Make check to verify that trailing edge is thick enough
[pn75, pn90] = low_thickness_exceedance(parameters, experiment_result);
if or(pn75 > 0, pn90 > 0 )
    % Limit pn factors to 0 (shape crossing handled separately)
    pn75 = min(pn75, 1);   % 0 means no exceedance, 1 means full exceedance
    pn90 = min(pn90, 1);   % 0 means no exceedance, 1 means full exceedance
    % Damp out expected L over D
    sc = sc * (1 - pn75).^2 * (1 - pn90).^2;
end

% Make check to verify that sides do not cross
valid_shape = valid_shape_test(experiment_result);
if not(valid_shape)
   sc = 0; 
end


% Make additional check to verify that L/D did not enter into unreasonable
% branch