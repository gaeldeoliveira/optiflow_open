%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Multiobjective Optimization Example Case for Airfoils with
%           Flaps. 
%           Arbitrary number of polar-configurations per design point.
%
%           Aerodynamic Goal: Glide ration L/D over a range of angles of 
%                             attack and flap angles
%           Structural  Goal: Maximum Airfoil Building Height
%
%           Parallel Execution Enabled
%
%           Resorts to include file: C1_include_FLAP_anysim_make_SD_SDD
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
    SC.N_cores = 4;                     % If parallelized, over how many cores (not cluster compatible)
%  ---- End of property edition ---- %
SC.set_context;

 %% Create Simulation Objects (Protocol and Workers)
 % Make a Simulation Protocol Object
SP = simulation_protocol('forced_transition' , SC);
%  ---- Edit any properties of choice ---- %
    SP.target_application='RBINKOLIVEIRA_V2';                      % To use Xfoil as your application set this to 'xfoil'
    SP.operation = 'alfa_polar_ref_start';                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
    SP.operation_parameters = [-15 20 0.5 0];                      % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
    SP.Re = 3000000;                                               % Reynolds Number (NaN means inviscid)
    SP.xtr_top=.99;                                                % set free transition on top
    SP.xtr_bot=.99;                                                % set free transition on bottom
    SP.N_crit=9;
%  ---- End of property edition ---- %


%% Create Shape and Parametrization handling objects
% Create Parametrization objects for upper and lower side
p_upper = parametrization( 'cst_upper'  , 8);
p_lower = parametrization( 'cst_lower'  , 8);
 
% Fix number of non-CST-shape parameters
N_dummy_parameters = 3;                                            % Number of parameters that are not shape definition ones for the three flap cases (undeflected, downward and upward)
N_flap_steps_per_side = 1;                                         % Number of flap steps per side 
 
%% Include Flap Simulation handling objects
C1_include_FLAP_anysim_make_SD_SDD

%% Declare Cost Functions (Single Objectives)
% Create a simple aerodynamic cost function
CFaero = cost_function('aerodynamic_sample');                      % Initialize it and give it a name
CFaero.cost_function_handle = @cost_function_FLAP_interp;          % Specify which function to use
CFaero.post_function = @(x) -x;                                     % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFaero.parameter_list = {'flap_angle_scalings'};                   % list of parameters names that control cost function execution. In this case we use it to pass information to distinguish 
CFaero.parameter_values  = {flap_angle_scalings};                  % list of parameters values that control cost function execution

% Create a simple structural cost function
CFstruct = cost_function('building height');                       % Initialize it and give it a name
CFstruct.cost_function_handle = @cost_function_thickness;          % Specify which function to use
CFstruct.post_function = @(x) -x;                                  % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFstruct.parameter_list = {};                                      % list of parameters names that control cost function execution
CFstruct.parameter_values  = {};                                   % list of parameters values that control cost function execution


%% Create the global cost function using the single objectives we defined earlier

% Instantiate a global_cost_function object
CFG = global_cost_function('global_cost_function_FLAP');                         
% Provide list of simulation workers to be executed on each cost function call    
CFG.simulation_worker_list = SW_cell_vector;
% Provide list of cost functions for interpretation of results
CFG.interpretation_function_list = {CFaero , CFstruct };
% Specify which experiment results (from SW) are passed to each interpretation function
CFG.argument_connectivity_matrix = [transpose(1:length(SW_cell_vector)) , ones(length(SW_cell_vector) , 1)];           
 
%% Set meaningful constraints on the design space
% Make Shape Fit Object for constraint suggestion , using previously defined Shape Definition object for the undeflected flap
SF = shape_fit_cst(SD, 'cst88_undeflected_flap');                  % 'fitcst88_undeflected_flap' is any arbitrary name you like

% Make a constraint manager object
CM = constraint_manager_cst('Constraint manager 1', SC, SD, SF);   % Create Constraint Manager Object, using SD1 (undeflected) shape definition as reference
%  ---- Edit any properties of choice ---- %
    % For example modify the extension factors that quantify how far we can get
    % from the smallest elementary einterval containing the reference designs
    CM.extension_factor_ub_LE = 1.8;                                % Greater value permits thicker airfoils
    CM.extension_factor_ub_M  = 1.5;                                % Greater value permits thicker airfoils
    CM.extension_factor_ub_TE = 2.5;                                % Greater value permits thicker airfoils
        
    CM.extension_factor_lb_LE = 1.4;                                % Greater value permits thinner airfoils
    CM.extension_factor_lb_M = 1.4;                                 % Greater value permits thinner airfoils
    CM.extension_factor_lb_TE = 1.4;                                % Greater value permits thinner airfoils

    % Use this to modify guessed values
    % CM.A_ineq = (...)

    % Check whether or not there are constraints on suction panel location
    CM.free_suction_area_constraints = false;                       % Do NOT Activate Constraints for free suction area
%  ---- End of property edition ---- %
CM.suggest_constraints();                                           % Update constraints based on option changes

%% Edit Flap Constraints Manually and Adapt input population

% Set Constraints Manually
CM.lb_ext(end-1) = 0.6 ; CM.ub_ext(end-1) = 0.8;                    % Set minimum and maximum hinge positions
CM.lb_ext(end) = 2 ; CM.ub_ext(end) = 10;                           % Set minimum and maximum flapping angle ranges

% Now make a random seeding of constraint-complying initial-population points for non-shape parameters
% For hinge positions first
CM.ref_pop_array_bounds(:, end-1) = rand(size(CM.ref_pop_array_bounds(:, end-1))) * (CM.ub_ext(end-1) - CM.lb_ext(end-1)) + CM.lb_ext(end-1);
% For flap ranges second
CM.ref_pop_array_bounds(:, end)   = rand(size(CM.ref_pop_array_bounds(:, end-1))) * (CM.ub_ext(end) - CM.lb_ext(end))     + CM.lb_ext(end);

%% Make genetic optimization manager object
GM = gamultiobj_manager('gamultiobj', CM, CFG);                     % Initialize it and tell it which optimizer to use (ga or gamultiobj)
%  ---- Edit any properties of choice ---- %
% Set the options controlling use of constraints and optimizer options
    GM.vectorize = false;                                           % Go Multicore! Vectorize == true means parallel execution becomes possible if system_context and simulation_workers are set for that. Serial execution is slower but better for preview. It is a good practice to set your case up as serial, and once activate vectorization once you know it works well.
    GM.use_reference_population = true;                             % Have a nice start!
%  ---- End of property edition ---- %
GM.set_options();                                                   % Create option structure according to above defined options
GM.options_structure.PopulationSize=70;
GM.options_structure.Generations=150;
GM.build_problem_structure();                                       % Create problem structure using above options and constraints from CM object
GM.options_structure.MutationFcn= gaoptimset(GM.options_structure, 'MutationFcn', {@mutationadaptfeasible, 0.05});

%% Run Optimization
% If you want to run your optimization without the GUI, uncomment the three
% following lines
GM.start_genetic_optimizer();
results = GM.last_results;
save(['results-', datestr(now, 'dd-mm-yyyy')], 'results');
 
 
 
 
 