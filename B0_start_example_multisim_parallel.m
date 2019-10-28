%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Multiobjective Optimization Example Case with Two Polars per
%           Cost Function Evaluation to enable explicit comparison of Free 
%           and Forced Transition performance

%           Aerodynamic Goal: Max Cl difference between Free and Force
%                             Transition Cases
%           Structural  Goal: Maximum Airfoil Thickness 
%
%           Parallel Execution Enabled
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all
 %clear all
 
%% System Work necessary to start

% Start by adding the source folders to the matlab path
fs = filesep();      % Folder separator is OS dependent
addpath([cd fs 'src']);
addpath([cd fs 'user_src']);
addpath([cd fs 'gui']);


% Create System Context Object and Set Context
SC = system_context;
%  ---- Edit any properties of choice, like number of cores ---- %
    SC.N_cores = 4;
    
%  ---- End of property edition ---- %
SC.set_context;


%% Create Shape and Parametrization handling objects
 % Create Parametrization objects for upper and lower side
 p_upper = parametrization( 'cst_upper'  , 8);
 p_lower = parametrization( 'cst_lower'  , 8);
 
 
 % Make Shape Definition Object
 % Create Shape Definition Objects using previously defined parametrizations
 N_dummy_parameters = 1;                                                        % Number of parameters that are not shape definition ones
 SD = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'cst44');    % 'cst44' is any arbitrary name you like
 
 % Make Shape Fit Object for constraint suggestion , using previously
 % defined Shape Definition object
 SF = shape_fit_cst(SD, 'fitcst88');                            % 'fitcst88' is any arbitrary name you like

%% Create Simulation Objects (Protocol and Worker)
 % Make a Simulation Protocol Object
SP1 = simulation_protocol('forced_transition' , SC);
%  ---- Edit any properties of choice ---- %
    SP1.target_application='RBINKOLIVEIRA_V2';                      % To use Xfoil as your application set this to 'xfoil'
    SP1.operation = 'alfa_polar_ref_start';                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
    SP1.operation_parameters = [-15 20 0.5 0];                      % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
    SP1.Re = 3000000;                                               % Reynolds Number (NaN means inviscid)
    SP1.xtr_top=.05;                                                % set forced transition on top
    SP1.xtr_bot=.10;                                                % set forced transition on bottom
    SP1.N_crit=9;
%  ---- End of property edition ---- %

 % Make a simulation worker for the previously defined Simulation Object
SW_1 = simulation_worker('simul_worker_1', SP1, [] , SD,  SC);   % Mariline
%  ---- Edit any properties of choice ---- %
    SW_1.app_name = 'RBINKOLIVEIRA_V2';                             % For xfoil write SW_1.app_name = 'xfoil'
    SW_1.parallelized=1;                                            % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
%  ---- End of property edition ---- %

SP2 = simulation_protocol('clean_transition' , SC);
%  ---- Edit any properties of choice ---- %
    SP2.target_application='RBINKOLIVEIRA_V2';                      % To use Xfoil as your application set this to 'xfoil'
    SP2.operation = 'alfa_polar_ref_start';                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
    SP2.operation_parameters = [-15 20 0.5 0];                      % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
    SP2.Re = 6000000;                                               % Reynolds Number (NaN means inviscid)
    SP2.xtr_top=.99;                                                % set forced transition on top
    SP2.xtr_bot=.99;                                                % set forced transition on bottom
    SP2.N_crit=9;
%  ---- End of property edition ---- %

 % Make a simulation worker for the previously defined Simulation Object
SW_2 = simulation_worker('simul_worker_2', SP2, [] , SD,  SC);   % Mariline
%  ---- Edit any properties of choice ---- %
    SW_2.app_name = 'RBINKOLIVEIRA_V2';                             % For xfoil write SW_1.app_name = 'xfoil'
    SW_2.parallelized=1;                                            % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
%  ---- End of property edition ---- %


%% Create Cost Functions (Single Objectives)
% Create a simple aerodynamic cost function
CFaero = cost_function('aerodynamic_sample');                       % Initialize it and give it a name
CFaero.cost_function_handle = @cost_function_Delta_Cl_Max;           % Specify which function to use
CFaero.post_function = @(x) abs(x);                             % post_function, applied to result of cost_function_handle, for example to invert values or make negative
% CFaero.parameter_list = {'cl_i' , 'w_i'};                        % list of parameters names that control cost function execution
% CFaero.parameter_values  = {0.5:0.1:1.2 , [0.05 0.1, 0.1 , 0.1 , 0.1 , 0.15 , 0.2 , 0.2 ]};           % list of parameters values that control cost function execution

% Create a simple structural cost function
CFstruct = cost_function('building height');                        % Initialize it and give it a name
CFstruct.cost_function_handle = @cost_function_thickness;           % Specify which function to use
CFstruct.post_function = @(x) -x;                                   % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFstruct.parameter_list = {};                                       % list of parameters names that control cost function execution
CFstruct.parameter_values  = {};                                    % list of parameters values that control cost function execution


%% Create the global cost function using the single objectives we defined earlier
% Use the global_cost_function class
CFG = global_cost_function('my_first_one');                     % Initialize it and give it a name
    CFG.simulation_worker_list = {SW_1, SW_2};                            % Provide list of simulation workers to be executed on each cost function call
    CFG.interpretation_function_list = {CFaero , CFstruct };        % Provide list of cost functions for interpretation of results
    CFG.argument_connectivity_matrix = [1 1 ; 2 1 ];                       % Specify how each cost function uses results from experiments
 
% CM.lb_ext = (...)
% CM.lb_0   = (...)

%% Set meaningful constraints on the design space
% Make a constraint manager object
CM = constraint_manager_cst('Constraint manager 1', SC, SD, SF);    % Create Constraint Manager Object
%  ---- Edit any properties of choice ---- %
    % For example modify the extension factors that quantify how far we can get
    % from the smallest elementary einterval containing the reference designs
    CM.extension_factor_ub_LE = 1.8;                            % Greater value permits thicker airfoils
    CM.extension_factor_ub_M  = 1.5;                            % Greater value permits thicker airfoils
    CM.extension_factor_ub_TE = 2.5;                            % Greater value permits thicker airfoils
        
    CM.extension_factor_lb_LE = 1.4;                            % Greater value permits thinner airfoils
    CM.extension_factor_lb_M = 1.4;                             % Greater value permits thinner airfoils
    CM.extension_factor_lb_TE = 1.4;                            % Greater value permits thinner airfoils

    % Use this to modify guessed values
    % CM.A_ineq = (...)

    % Check whether or not there are constraints on suction panel location
    CM.free_suction_area_constraints = false;   % Do NOT Activate Constraints for free suction area
%  ---- End of property edition ---- %
CM.suggest_constraints();                                       % Update constraints based on option changes

%% Make genetic optimization manager object
GM = gamultiobj_manager('gamultiobj', CM, CFG);                 % Initialize it and tell it which optimizer to use (ga or gamultiobj)
%  ---- Edit any properties of choice ---- %
% Set the options controlling use of constraints and optimizer options
    GM.vectorize                    = true;                     % Slower but better for preview. You should always set your case up this way, and once you know it works well, activate vectorization.
    GM.use_reference_population     = true;                     % Have a nice start!
%  ---- End of property edition ---- %
GM.set_options();                                               % Create option structure according to above defined options
GM.options_structure.PopulationSize = 96;                       % Small population size for faster setup (2 * 3 batches of 32 cores)
GM.options_structure.Generations    = 50;                       % Small generation count for first results
GM.build_problem_structure();            % Create problem structure using above options and constraints from CM object
GM.options_structure.MutationFcn= gaoptimset(GM.options_structure, 'MutationFcn', {@mutationadaptfeasible, 0.05});
 
%% Run Optimization
% If you want to run your optimization without the GUI, uncomment the three
% following lines
GM.start_genetic_optimizer();
results = GM.last_results;
save(['results-', datestr(now, 'dd-mm-yyyy')], 'results');
 
 
 
 
 
 