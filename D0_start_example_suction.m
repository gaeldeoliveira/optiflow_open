%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Basic Multiobjective Optimization Example Case using Boundary
%           Layer Suction.
%           Simulation Protocol Dynamizer Object allows arbitrary suction
%           panel position and lenght, with fixed Cq.
%
%           Aerodynamic Goal: L/D over a range of Cls
%           Structural  Goal: Building Height
%
%           Parallel Execution Enabled
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all ; clear all ; clc ;

%% System Work necessary to start

% Start by adding the source folders to the matlab path
if ispc             % Folder separator is OS dependent
    fs = '\';       % For windows
else
    fs = '/';       % For mac or linux
end
addpath([cd fs 'src']);
addpath([cd fs 'user_src']);
addpath([cd fs 'gui']);


% Create System Context Object and Set Context
SC = system_context;
%  ---- Edit any properties of choice, like number of cores ---- %
    SC.N_cores = 1;
    
    % Create a temporary directory based on either handles.casename for GUI
    % launch casename
    if exist('handles'); SC.tmp_subdir = [ 'tmp' handles.casename fs]; end %#ok<EXIST>
    % Or based on time otherwise
    if not(exist('handles')); SC.tmp_subdir = [ 'tmp' datestr(now, 30) fs]; end %#ok<EXIST>
    
%  ---- End of property edition ---- %
SC.set_context;

%% Shape Definition
% Create Parametrization objects for upper and lower side
p_upper = parametrization( 'cst_upper'  , 8);
p_lower = parametrization( 'cst_lower'  ,  8);

% Make Shape Definition Object
% Create Shape Definition Objects using previously defined parametrizations
N_dummy_parameters = 2;
SD = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'cst66');

% Make Shape Fit Object for constraint suggestion , using previously
% defined Shape Definition object
SF = shape_fit_cst(SD, 'fitcst66');


%% Aerodynamic Simulation Protocol
% Make a Simulation Object
SP1 = simulation_protocol('clean_no_suction' , SC);
%  ---- Edit any properties of choice ---- %
SP1.target_application = 'rfoilsuc';                   

SP1.operation = 'alfa_polar_ref_start';
SP1.operation_parameters = [-10 20 0.5 0];        % Make polar for 0 to 20 degrees in 1 degree increments

SP1.Re = 3e6;             % Reynolds Number (NaN means inviscid)

% Simulate roughness. Forced transition on upper side at 1% and lower side at 5%
SP1.xtr_top = 0.01;           % x/c of forced transition point for top side (1 = free)
SP1.xtr_bot = 0.05;           % x/c of forced transition point for bottom side (1 = free)

SP1.suction_distribution_parameters.ssuc = 1;               % Top Side
SP1.suction_distribution_parameters.xsuc = [0.7 0.9];       % Suction panel from 70 to  90%
SP1.suction_distribution_parameters.vsuc = -0.003;          % Suction speed (to Uinf)
SP1.suction_distribution_type = 'rfoil_NSM2';               % Data protocol for suction interfacing
%  ---- End of property edition ---- %

% Make a Simulation Dynamizer for the previewsly defined Simulation Object
SPD1 = simulation_protocol_dynamizer('Fixed_Suction_Mass', SP1);               % This is used to handle the case in which the optimizer also searches for the best suction region
SPD1.fixed_Cqm = -0.004*0.2;
SPD1.ssuc = 1;

% Make a simulation worker for the previously defined Simulation Object
SW_1 = simulation_worker('simul_worker_1', SP1, SPD1 , SD,  SC);   % Mariline
SW_1.app_name = 'RBINKOLIVEIRA_V2';

% Create a simple aerodynamic cost function
CFaero = cost_function('aerodynamic_sample');                       % Initialize it and give it a name
CFaero.cost_function_handle = @cost_function_CLCD_at_CL;           % Specify which function to use
CFaero.post_function = @(x) abs(1 / x);                             % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFaero.parameter_list = {'cl_i' , 'w_i'};                        % list of parameters names that control cost function execution
CFaero.parameter_values  = {0.5:0.1:1.2 , [0.05 0.1, 0.1 , 0.1 , 0.1 , 0.15 , 0.2 , 0.2 ]};           % list of parameters values that control cost function execution

% Create a simple structural cost function
CFstruct = cost_function('building height');                        % Initialize it and give it a name
CFstruct.cost_function_handle = @cost_function_thickness;           % Specify which function to use
CFstruct.post_function = @(x) -x;                                   % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFstruct.parameter_list = {};                                       % list of parameters names that control cost function execution
CFstruct.parameter_values  = {};                                    % list of parameters values that control cost function execution


% Now make global cost function using previously defined single objective
% cost functions
CFG = global_cost_function('my_first_one');                         % Initialize it and give it a name
CFG.simulation_worker_list = {SW_1};                                % Provide list of simulation workers to be executed on each cost function call
CFG.interpretation_function_list = {CFaero , CFstruct};             % Provide list of cost functions for interpretation of results
CFG.argument_connectivity_matrix = [1 1];                           % Specify how each cost function uses results from experiments
 
% % Make constraint manager object
CM = constraint_manager_cst('Constraint amanger 1', SC, SD, SF);    % Create Constraint Manager Object
%  ---- Edit any properties of choice ---- %
% Use this to modify guessed values
% CM.A_ineq = (...)
% CM.lb_ext = (...)
% CM.lb_0   = (...)

CM.extension_factor_ub_LE = 1.5;            % Greater value permits thicker airfoils
CM.extension_factor_ub_M  = 1.3;            % Greater value permits thicker airfoils
CM.extension_factor_ub_TE = 2;            % Greater value permits thicker airfoils
        
CM.extension_factor_lb_LE = 1.2;            % Greater value permits thinner airfoils
CM.extension_factor_lb_M = 1.2;             % Greater value permits thinner airfoils
CM.extension_factor_lb_TE = 1.2;            % Greater value permits thinner airfoils

CM.x_min_suction = 0.4;                                             % Suction can start at 50% chord
CM.free_suction_area_constraints = true;                            % Activate Constraints for free suction area
%  ---- End of property edition ---- %
CM.suggest_constraints();                                           % Update constraints based on option changes

% Make genetic optimization manager object
GM = gamultiobj_manager('gamultiobj', CM, CFG);                     % Initialize it and tell it which optimizer to use (ga or gamultiobj)
%  ---- Edit any properties of choice ---- %
% Set the options controlling use of constraints and optimizer options
GM.vectorize = false;                                               % Slower but better for preview. You should always set your case up this way, and once you know it works well, activate vectorization.
%  ---- End of property edition ---- %
GM.set_options();                                                   % Create option structure according to above defined options
GM.build_problem_structure();                                       % Create problem structure using above options and constraints from CM object

%% Run Optimization
% If you want to run your optimization without the GUI, uncomment the three
% following lines
GM.start_genetic_optimizer();
results = GM.last_results;
save(['results-', datestr(now, 'dd-mm-yyyy')], 'results');

%% Obsolete code snippets left for later reference!
% % Experience with linear constraints
% CM.A_ineq = zeros(2 , 21);
% CM.A_ineq(1:2 , 20:21) = [[1 -1] ; [-1 1]];
% CM.b_ineq = zeros(2, 1);
% CM.b_ineq(1) = 0.1;         % Minimum Length of suction area
% CM.b_ineq(2) = 0.3;         % Maximum Length of suction area
% 
% GM.use_ineq_constraints = false;
% 
% GM.set_options();                                                   % Create option structure according to above defined options
% GM.build_problem_structure();                                       % Create problem structure using above options and constraints from CM object

% GM.start_genetic_optimizer();