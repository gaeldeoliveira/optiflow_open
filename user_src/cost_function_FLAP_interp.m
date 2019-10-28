function sc = cost_function_FLAP_interp(parameters , experiments_results)
%sc = cost_function_FLAP_interp(parameters , experiments_results)
%   Demonstrates the construction of a 2d scattered interpolant to construct 
%       alpha-beta polar curves for airfoils with control surfaces (such as flaps!)!
%       A similar structure could also be used for multi-point optimization
%       of other actuation mechanisms, for example based on the multi 

%% Extract the data we need to build interpolants
% Extract parameters from parameters structure
    % No parameters for this case

% Extract Relevant parameters from static parameters structure
flap_angle_scalings = parameters.flap_angle_scalings; 

% Identify number of experiments and flapping angle range
N_experiments = length(experiments_results);
beta_range = experiments_results{1}.x(end);     % x is the same for all experiments/samples in a set!

% Declare scattered interpolant point vectors
alpha_vector =[];
beta_vector  =[];
cl_vector    =[];
cd_vector    =[];
cm_vector    =[];

% Run over complete set of experiments to extract data to construct objects 
% from standard input form
for n_experiment = 1:N_experiments
    % First check there exists a converged and properly loaded polar for
    % each experiment
    % Determine beta range for this experiment
    
    if strcmp(class(experiments_results{n_experiment}.ap) , 'aerodynamic_polar')
        % Extract alpha, Cl, Cm, Cd
        alpha_experiment   = experiments_results{n_experiment}.ap.alpha_range();
        cl_experiment      = experiments_results{n_experiment}.ap.cl_alpha(alpha_experiment);        
        cd_experiment      = experiments_results{n_experiment}.ap.cd_alpha(alpha_experiment);
        cm_experiment      = experiments_results{n_experiment}.ap.cm_alpha(alpha_experiment);
        
        % Store into interpolant vectors
        alpha_vector       = [alpha_vector  ; alpha_experiment(:)];             %#ok<AGROW>
        cl_vector          = [cl_vector     ; cl_experiment(:)];                %#ok<AGROW>
        cd_vector          = [cd_vector     ; cd_experiment(:)];                %#ok<AGROW>
        cm_vector          = [cm_vector     ; cm_experiment(:)];                %#ok<AGROW>
        
        % Generate beta interpolant vector for this experiment and append
        % to the large one
        % Deprecate:
        %       beta_experiment = ones(size(alpha_experiment(:))) * flap_angle_scalings(n_experiment) * beta_range;
        % And switch to cleaner approach
        beta_experiment = ones(size(alpha_experiment(:))) * experiments_results{n_experiment}.multisim_id_value * beta_range;
        beta_vector  =[beta_vector ; beta_experiment];                          %#ok<AGROW>
        
        
        % Proceed to next experiment!
    end 
end

%% Construct the interpolants to get robust alpha-beta polars for cl, cd and cm

interpolation_method = 'linear';              % Linear offers garanteed robustness, but 'natural' could also be interesting as it is C1 continuous except at data points  
extrapolation_method = 'none';                % Return NaN whenever we are outside of the convex hull and deal with it afterwards
% First make a classical scattered interpolant returning NaN's 
clInterpolant = scatteredInterpolant(alpha_vector(:), beta_vector(:), cl_vector(:), interpolation_method , extrapolation_method);
cdInterpolant = scatteredInterpolant(alpha_vector(:), beta_vector(:), cd_vector(:), interpolation_method , extrapolation_method);
cmInterpolant = scatteredInterpolant(alpha_vector(:), beta_vector(:), cm_vector(:), interpolation_method , extrapolation_method);

% And now complement ou scattered interpolant with a NaN handling
cl_fallback_value = 0;      % Verify these values always introduce a cost function value degradation compared to converged elements
cd_fallback_value = 1;
cm_fallback_value = 0;      % Cm needs to be checked!

% Make interpolant functions using the fallbacks (a=alpha, b=beta for shorter line!)
cl_alpha_beta = @(a,b) nan_to_zero(clInterpolant(a,b)) + isnan(clInterpolant(a,b)) .* cl_fallback_value;
cd_alpha_beta = @(a,b) nan_to_zero(cdInterpolant(a,b)) + isnan(cdInterpolant(a,b)) .* cd_fallback_value;
cm_alpha_beta = @(a,b) nan_to_zero(cmInterpolant(a,b)) + isnan(cmInterpolant(a,b)) .* cm_fallback_value;


%% Now compute our cost function (this is a very very simplistic approach, and in fact the interpolant let's you integrate over any AOA/Flap history you may want!
AOA_mean    =  5; 
AOA_range   = 15;

sc = cl_alpha_beta(AOA_mean+AOA_range/2, +beta_range) - cl_alpha_beta(AOA_mean-AOA_range/2, -beta_range);


%% 
if or(isnan(sc), isinf(sc))
    % If cost function got NaN or Inf by some strange behavior anywhere throw it out !
    sc = 0;
end