%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Basic Multiobjective Optimization Example Case
%
%           Aero Goal 1: max(expected(L)/expected(D)) with free transition
%           Aero Goal 2: max(expected(L)/expected(D)) with trip transition
%
%           Structural constraint: t/c = 0.24 (Inwind tip)
%
%           Parallel Execution Enabled
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all; clc;
%
% Based on caseB_Re9_alpha7_tc21_sigma05

%
%
%
%
 
%% Case definition
case_overview.N_upper            =   8   ;       % [int. ] - order of upper side parametrization (order=degree+1)
case_overview.N_lower            =   8   ;       % [int. ] - order of lower side parametrization (order=degree+1)
case_overview.N_dummy_parameters =   0   ;       % [int. ] - number of dummy (non-shape) parameters in genotype vector
case_overview.target_thickness   =   0.24;       % [adim.] - thickness over chord of optimized airfoils (reprojection)

case_overview.polar_range    = [-15 22 0.2 0];   % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
case_overview.Re             = 9000000   ;       % Reynolds Number (NaN means inviscid)
case_overview.xtr_top_free       =   0.99;       % set top    transition position for free    case
case_overview.xtr_bot_free       =   0.99;       % set bottom transition position for free    case
case_overview.xtr_top_trip       =   0.05;       % set top    transition position for tripped case
case_overview.xtr_bot_trip       =   0.10;       % set bottom transition position for tripped case
case_overview.N_crit             =   9   ;       % critical amplification factor for Tollmien-Schlichting waves

case_overview.N_cores            =   2   ;       % [int. ] - number of cores over which to spread computations
case_overview.parallelized       =   1   ;       % [ibool] - parallelize simulations of cost function computation (requires that vectorize is true to be effective)
case_overview.vectorized         =  true ;       % [bool.] - allow genetic algorithm to make vectorized cost function calls

case_overview.PopulationSize     =  90   ;       % [int. ] - size of population in genetic pool
case_overview.Generations        =  30   ;       % [int. ] - number of generations to run genetic algorithm

case_overview.alpha_bar          =   7   ;       % [deg. ] - mean angle of attack (around which expected values are computed)              
case_overview.n_alpha_til        =  64   ;       % [int. ] - number of elements convolution integral kernel
case_overview.sigma              =   0.5 ;       % [deg. ] - standard deviation of angle of attack fluctuations
case_overview.n_sigma            =   1.64485;    % [adim.] - Width of integration range (in sigmas, affects pdf/cdf coverage)
case_overview.n_free_exp         =   1   ;       % [int. ] - Number of experiment (or simulation worker-protocol pair) corresponding to free   transition polar
case_overview.n_trip_exp         =   2   ;       % [int. ] - Number of experiment (or simulation worker-protocol pair) corresponding to forced transition polar

%% System Work necessary to start
% Start by adding the source folders to the matlab path
fs = filesep();                                                                     % Folder separator is OS dependent
addpath([cd fs 'src']);                                                             % Add optimizer source code folder
addpath([cd fs 'user_src']);                                                        % Add optimizer user sources (cost functions, etc) folder
addpath([cd fs 'gui']);                                                             % Add gui source code folder
addpath([cd fs 'dev']);                                                             % Add gui source code folder
addpath([cd fs 'probabilistic']);                                                   % Add gui source code folder


% Create System Context Object and Set Context
SC = system_context;
%  ---- Edit any properties of choice, like number of cores ---- %
    SC.N_cores = case_overview.N_cores;                                             % N_cores is actually equivalent to number of Threads. (this was written when cores were single threaded, on a Core 2 Duo ULV!)
    SC.airfoil_subdir             = 'probabilistic/data/airfoils/';
    SC.feasible_airfoil_subdir = 'probabilistic/data/feasible_airfoils/';
%  ---- End of property edition ---- %
SC.set_context;


%% Create Shape and Parametrization handling objects
 % Create Parametrization objects for upper and lower side
 p_upper = parametrization( 'cst_upper'  , case_overview.N_upper);
 p_lower = parametrization( 'cst_lower'  , case_overview.N_lower);
 
 % Make Shape Definition Object
 % Create Shape Definition Objects using previously defined parametrizations
 SD = shape_definition_cst(p_upper , p_lower , case_overview.N_dummy_parameters, 'cst44');    % 'cst44' is any arbitrary name you like
 
 % Make Shape Fit Object for constraint suggestion , using previously
 % defined Shape Definition object
 SF = shape_fit_cst(SD, 'fitcst88');                                                % 'fitcst88' is any arbitrary name you like

%% Create Simulation Objects (Protocol and Worker)
 % Make a Simulation Protocol Object
SP1 = simulation_protocol('free_transition' , SC);
%  ---- Edit any properties of choice ---- %
    SP1.target_application='RBINKOLIVEIRA_V2';                                      % To use Xfoil as your application set this to 'xfoil'
    SP1.operation = 'alfa_polar_ref_start';                                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
    SP1.operation_parameters = case_overview.polar_range;                           % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
    SP1.Re                   = case_overview.Re;                                    % Reynolds Number (NaN means inviscid)
    SP1.xtr_top              = case_overview.xtr_top_free;                          % set forced transition on top
    SP1.xtr_bot              = case_overview.xtr_bot_free;                          % set forced transition on bottom
    SP1.N_crit               = case_overview.N_crit;                                % critical amplification factor for Tollmien-Schlichting waves
%  ---- End of property edition ---- %

 % Make a simulation worker for the previously defined Simulation Object
SW_1 = simulation_worker('simul_worker_1', SP1, [] , SD,  SC);                      % Mariline
%  ---- Edit any properties of choice ---- %
    SW_1.app_name            = 'RBINKOLIVEIRA_V2';                                  % For xfoil write SW_1.app_name = 'xfoil'
    SW_1.parallelized        = case_overview.parallelized;                          % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
    SW_1.fig_active          = 1;                                                   % Activate plotting for clean case
%  ---- End of property edition ---- %

SP2 = simulation_protocol('tripped_transition' , SC);
%  ---- Edit any properties of choice ---- %
    SP2.target_application='RBINKOLIVEIRA_V2';                                      % To use Xfoil as your application set this to 'xfoil'
    SP2.operation = 'alfa_polar_ref_start';                                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
    SP2.operation_parameters = case_overview.polar_range;                           % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
    SP2.Re                   = case_overview.Re;                                    % Reynolds Number (NaN means inviscid)
    SP2.xtr_top              = case_overview.xtr_top_trip;                          % set forced transition on top
    SP2.xtr_bot              = case_overview.xtr_bot_trip;                          % set forced transition on bottom
    SP2.N_crit               = case_overview.N_crit;                                % critical amplification factor for Tollmien-Schlichting waves
%  ---- End of property edition ---- %

 % Make a simulation worker for the previously defined Simulation Object
SW_2 = simulation_worker('simul_worker_2', SP2, [] , SD,  SC);                      % Mariline
%  ---- Edit any properties of choice ---- %
    SW_2.app_name     = 'RBINKOLIVEIRA_V2';                                         % For xfoil write SW_1.app_name = 'xfoil'
    SW_2.parallelized = case_overview.parallelized;                                 % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
    SW_1.fig_active   = 0;                                                          % De-activate plotting for tripped case
%  ---- End of property edition ---- %

%% Create Cost Functions (Single Objectives)
% Define cost function parameters
alpha_bar   = case_overview.alpha_bar ;                                              % [deg. ] - mean angle of attack (around which expected values are computed)              
n_alpha_til = case_overview.n_alpha_til;                                             % [int. ] - number of elements convolution integral kernel
sigma       = case_overview.sigma;                                                   % [deg. ] - standard deviation of angle of attack fluctuations
n_sigma     = case_overview.n_sigma;                                                 % [adim.] - Width of integration range (in sigmas, affects pdf/cdf coverage)
n_free_exp  = case_overview.n_free_exp;                                              % [int. ] - Number of experiment (or simulation worker-protocol pair) corresponding to free   transition polar
n_trip_exp  = case_overview.n_trip_exp;                                              % [int. ] - Number of experiment (or simulation worker-protocol pair) corresponding to forced transition polar
ld_bar_ref_free  = 1;                                                                % [adim.] - Initial dummy value for L over D normalization 
cl_bar_ref_free  = 1;                                                                % [adim.] - Initial dummy value for Cl       normalization 
ld_bar_ref_trip  = 1;                                                                % [adim.] - Initial dummy value for L over D normalization 
cl_bar_ref_trip  = 1;                                                                % [adim.] - Initial dummy value for Cl       normalization 
bh_75_threshold  = 0;                                                                % [adim.] - Initial dummy value for minimum thickness threshold at 0.75c
bh_90_threshold  = 0;                                                                % [adim.] - Initial dummy value for minimum thickness threshold at 0.90c
% Create clean aerodynamics cost function
CFfree = cost_function('expected_free_LD_at_AoA');                                   % Initialize it and give it a name
CFfree.cost_function_handle = @cost_function_max_expected_CL_over_expected_CD;   % Specify which function to use
CFfree.post_function    = @(x) x;                                                    % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFfree.parameter_list   = {'alpha_bar', 'n_alpha_til', 'n_sigma', 'n_experiment', 'ld_bar_ref'    , 'cl_bar_ref'    , 'bh_75_threshold', 'bh_90_threshold', 'sigma'};   % list of parameters names that control cost function execution
CFfree.parameter_values = { alpha_bar ,  n_alpha_til ,  n_sigma ,  n_free_exp ,  ld_bar_ref_free ,  cl_bar_ref_free ,  bh_75_threshold ,  bh_90_threshold ,  sigma };   % list of parameters values that control cost function execution

% Create forced aerodynamics cost function
CFtrip = cost_function('expected_tripped_LD_at_AoA');                                % Initialize it and give it a name
CFtrip.cost_function_handle = @cost_function_max_expected_CL_over_expected_CD;   % Specify which function to use
CFtrip.post_function    = @(x) x;                                                    % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFtrip.parameter_list   = {'alpha_bar', 'n_alpha_til', 'n_sigma', 'n_experiment', 'ld_bar_ref'    , 'cl_bar_ref'    , 'bh_75_threshold', 'bh_90_threshold', 'sigma'};   % list of parameters names that control cost function execution
CFtrip.parameter_values = { alpha_bar ,  n_alpha_til ,  n_sigma ,  n_trip_exp ,  ld_bar_ref_trip ,  cl_bar_ref_trip ,  bh_75_threshold ,  bh_90_threshold ,  sigma };   % list of parameters values that control cost function execution


%% Create the global cost function using the single objectives we defined earlier
% Use the global_cost_function class
CFG = global_cost_function('Free_and_forced_expected_LD_at_AoA');                    % Initialize it and give it a name
CFG.simulation_worker_list = {SW_1, SW_2};                                           % Provide list of simulation workers to be executed on each cost function call
CFG.interpretation_function_list = {CFfree , CFtrip};                                % Provide list of cost functions for interpretation of results
CFG.argument_connectivity_matrix = [1 1 ; 2  2];                                     % Specify how each cost function uses results from experiments
 
% CM.lb_ext = (...)
% CM.lb_0   = (...)

%% Set meaningful constraints on the design space
% Make a constraint manager object
CM = constraint_manager_cst('Constraint manager 1', SC, SD, SF);                     % Create Constraint Manager Object
%  ---- Edit any properties of choice ---- %
    % For example modify the extension factors that quantify how far we can get
    % from the smallest elementary einterval containing the reference designs
    CM.extension_factor_ub_LE = 1.8;                                                 % Greater value permits thicker airfoils
    CM.extension_factor_ub_M  = 1.5;                                                 % Greater value permits thicker airfoils
    CM.extension_factor_ub_TE = 2.5;                                                 % Greater value permits thicker airfoils
        
    CM.extension_factor_lb_LE = 1.4;                                                 % Greater value permits thinner airfoils
    CM.extension_factor_lb_M  = 1.4;                                                 % Greater value permits thinner airfoils
    CM.extension_factor_lb_TE = 1.4;                                                 % Greater value permits thinner airfoils

    % Use this to modify guessed values
    % CM.A_ineq = (...)

    % Check whether or not there are constraints on suction panel location
    CM.free_suction_area_constraints = false;                                        % Do NOT Activate Constraints for free suction area
%  ---- End of property edition ---- %
CM.suggest_constraints();                                                            % Update constraints based on option changes

% %% Handle additional stuff for probabilistic design (if sigma is a design variable)
% % Adapt lower and upper bounds
% CM.lb_ext(end) = 0;
% CM.ub_ext(end) = 8.0;
% % Adapt reference population to include a constant turbulence intensity
% % CM.ref_pop_array_bounds(:,end) = 0.1;
% % Adapt reference population to include a randomized turbulence intensity
% a=CM.lb_ext(end); b=CM.ub_ext(end);                                                 % Get states back
% rng(25041974);                                                                      % Set random generator state with a good seed!
% r = a + (b-a).*rand(size(CM.ref_pop_array_bounds, 1),1);                            % Generate a uniformly sampled set of random numbers within bounds
% CM.ref_pop_array_bounds(:,end) = r;                                                 % Pass it into ref_pop_array_bounds (used for 


%% Make genetic optimization manager object
GM = gamultiobj_manager('gamultiobj', CM, CFG);                                       % Initialize it and tell it which optimizer to use (ga or gamultiobj)
%  ---- Edit any properties of choice ---- %
% Set the options controlling use of constraints and optimizer options
    GM.vectorize = case_overview.vectorized;
    GM.use_reference_population = true;                                               % Have a nice start!
    %GM.vectorize = false;                                                            % Slower but better for preview. You should always set your case up this way, and once you know it works well, activate vectorization.
    %GM.vectorize = true;                                                             % Go Multicore!
%  ---- End of property edition ---- %
GM.set_options();                                                                     % Create option structure according to above defined options
GM.options_structure.PopulationSize = case_overview.PopulationSize;
GM.options_structure.Generations    = case_overview.Generations   ;
GM.build_problem_structure();                                                         % Create problem structure using above options and constraints from CM object
GM.options_structure.MutationFcn= gaoptimset(GM.options_structure, 'MutationFcn', {@mutationadaptfeasible, 0.05});

 %% Make SDD Object for Thickness Projection (after constraint suggestion)
SDD = shape_dynamizer(SD);
SDD.mode = 'thickness';
SDD.target_thickness = case_overview.target_thickness;

%% Look at some reference airfoils
% Set starting airfoil (DU97W300, depends on state of data folders! WARNING)
%    7  AA21
%   10  DU21
%   12  FFA21
%   17  RISO21
%   20  S809

% Now select a small reference pool
n_reference                   = [7 10 12 17];
x_reference                   = GM.CM.ref_pop_array_bounds(n_reference ,:);
names_reference               = GM.CM.names_fitted_bounds(n_reference);
% And display its names
for n=1:length(names_reference); disp([num2str(n) ' ' names_reference{n}{1}]); end
% Compute cost function on it
cf_reference                  = GM.call_cost_function(x_reference);
% Now extract experiment results
reference_experiment_set      = GM.CFG.last_experiment_set;
% And corresponding parameters
reference_parameters_free     = CFfree.parameter_structure;
reference_parameters_trip     = CFtrip.parameter_structure;
% Allocate reference values for ld_bar and cl_bar
reference_ld_bar_free         = zeros(size(reference_experiment_set,1), 1);
reference_cl_bar_free         = zeros(size(reference_experiment_set,1), 1);
reference_ld_bar_trip         = zeros(size(reference_experiment_set,1), 1);
reference_cl_bar_trip         = zeros(size(reference_experiment_set,1), 1);
reference_bh75                = zeros(size(reference_experiment_set,1), 1);
reference_bh90                = zeros(size(reference_experiment_set,1), 1);
% Compute reference values of ld_bar and cl_bar for each airfoil
for n_ref = 1:length(reference_experiment_set)
    % Free
    reference_ld_bar_free(n_ref) = cost_function_expected_CLCD_at_expected_AoA(reference_parameters_free , reference_experiment_set(n_ref, :));
    reference_cl_bar_free(n_ref) = cost_function_expected_CL_at_expected_AoA(  reference_parameters_free , reference_experiment_set(n_ref, :));
    % Tripped
    reference_ld_bar_trip(n_ref) = cost_function_expected_CLCD_at_expected_AoA(reference_parameters_trip , reference_experiment_set(n_ref, :));
    reference_cl_bar_trip(n_ref) = cost_function_expected_CL_at_expected_AoA(  reference_parameters_trip , reference_experiment_set(n_ref, :));
    % Size
    [~, ~, reference_bh75(n_ref), reference_bh90(n_ref)] = low_thickness_exceedance(reference_parameters_trip , reference_experiment_set{n_ref,1});
end
% Now select best values for each 
ld_bar_ref_free = max(reference_ld_bar_free);
cl_bar_ref_free = max(reference_cl_bar_free);
% Tripped
ld_bar_ref_trip = max(reference_ld_bar_trip);
cl_bar_ref_trip = max(reference_cl_bar_trip);
% Set minimum thickness thresholds
bh_75_threshold = 0.60 * min(reference_bh75);
bh_90_threshold = 0.80 * min(reference_bh90);

% And fill them into the parameter structure (update)
CFfree.parameter_list   = {'alpha_bar', 'n_alpha_til', 'n_sigma', 'n_experiment', 'ld_bar_ref'    , 'cl_bar_ref'    , 'bh_75_threshold', 'bh_90_threshold', 'sigma'};   
CFfree.parameter_values = { alpha_bar ,  n_alpha_til ,  n_sigma ,  n_free_exp ,  ld_bar_ref_free ,  cl_bar_ref_free ,  bh_75_threshold ,  bh_90_threshold ,  sigma };
CFtrip.parameter_list   = {'alpha_bar', 'n_alpha_til', 'n_sigma', 'n_experiment', 'ld_bar_ref'    , 'cl_bar_ref'    , 'bh_75_threshold', 'bh_90_threshold', 'sigma'};   
CFtrip.parameter_values = { alpha_bar ,  n_alpha_til ,  n_sigma ,  n_trip_exp ,  ld_bar_ref_trip ,  cl_bar_ref_trip ,  bh_75_threshold ,  bh_90_threshold ,  sigma };   % list of parameters values that control cost function execution



% Test cost function anew!
cf_new_reference        = GM.call_cost_function(x_reference);


% %% Now, start optimizer (with phenotype space skewing)
% disp('-- Ready to Start Optimization --')
% GM.skew_phenotype_space        = true;

%% Test cost function
GM.call_cost_function(GM.CM.ref_pop_array_bounds(2,:));
 
%% Run Optimization
% To run optimization without the GUI, uncomment the four following lines
% GM.start_genetic_optimizer();
% results = GM.last_results;
% save(['probabilistic/results/', mfilename, '_results_' , datestr(now, 'yyyymmddTHHMMSS')], 'results');
% save(['probabilistic/results/', mfilename, '_fullcase_', datestr(now, 'yyyymmddTHHMMSS')]);

    



