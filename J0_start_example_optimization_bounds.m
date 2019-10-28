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
%           Aerodynamic Goal: L/D over a range of Cls
%           Structural  Goal: Shell stiffness (Rody Kemp, under Gael)
%
%           Parallel Execution Enabled
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all
 %clear all
 
%% System Work necessary to start
% Start by adding the source folders to the matlab path
fs = filesep();                 % Folder separator is OS dependent
addpath([cd fs 'src']);         % Add optimizer source code folder
addpath([cd fs 'user_src']);    % Add optimizer user sources (cost functions, etc) folder
addpath([cd fs 'gui']);         % Add gui source code folder


% Create System Context Object and Set Context
SC = system_context;
%  ---- Edit any properties of choice, like number of cores ---- %
    SC.N_cores = 1;             % N_cores is actually equivalent to number of Threads. (this was written when cores were single threaded on a Core 2 Duo ULV!)
    
%  ---- End of property edition ---- %
SC.set_context;


%% Create Shape and Parametrization handling objects
 % Create Parametrization objects for upper and lower side
 p_upper = parametrization( 'cst_upper'  , 6);
 p_lower = parametrization( 'cst_lower'  , 6);
 
 
 % Make Shape Definition Object
 % Create Shape Definition Objects using previously defined parametrizations
 N_dummy_parameters = 0;                                                        % Number of parameters that are not shape definition ones
 SD = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'cs66');    % 'cst44' is any arbitrary name you like
 
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
    SP1.Re = 6000000;                                               % Reynolds Number (NaN means inviscid)
    SP1.xtr_top=.99;                                                % set forced transition on top
    SP1.xtr_bot=.99;                                                % set forced transition on bottom
    SP1.N_crit=9;
%  ---- End of property edition ---- %

 % Make a simulation worker for the previously defined Simulation Object
SW_1 = simulation_worker('simul_worker_1', SP1, [] , SD,  SC);   % Mariline
%  ---- Edit any properties of choice ---- %
    SW_1.app_name = 'RBINKOLIVEIRA_V2';                             % For xfoil write SW_1.app_name = 'xfoil'
    SW_1.parallelized=1;                                            % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
    SW_1.SC.N_cores=4;                                              % If parallelized, over how many cores (not cluster compatible)
%  ---- End of property edition ---- %

%% Create Cost Functions (Single Objectives)
% Create a simple aerodynamic cost function
CFaero = cost_function('aerodynamic_sample');                       % Initialize it and give it a name
CFaero.cost_function_handle = @cost_function_CLCD_at_CL;           % Specify which function to use
CFaero.post_function = @(x) abs(1 / x);                             % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFaero.parameter_list = {'cl_i' , 'w_i'};                        % list of parameters names that control cost function execution
CFaero.parameter_values  = {0.5:0.1:1.2 , [0.05 0.1, 0.1 , 0.1 , 0.1 , 0.15 , 0.2 , 0.2 ]};           % list of parameters values that control cost function execution

% Create a simple structural cost function
CFstruct = cost_function('building height');                        % Initialize it and give it a name
CFstruct.cost_function_handle = @cf_free_StifnessFinal;           % Specify which function to use
CFstruct.post_function = @(x) -x;                                   % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFstruct.parameter_list = {};                                       % list of parameters names that control cost function execution
CFstruct.parameter_values  = {};                                    % list of parameters values that control cost function execution


%% Create the global cost function using the single objectives we defined earlier
% Use the global_cost_function class
CFG = global_cost_function('my_first_one');                     % Initialize it and give it a name
    CFG.simulation_worker_list = {SW_1};                            % Provide list of simulation workers to be executed on each cost function call
    CFG.interpretation_function_list = {CFaero , CFstruct };        % Provide list of cost functions for interpretation of results
    CFG.argument_connectivity_matrix = [1 1];                       % Specify how each cost function uses results from experiments
 
% CM.lb_ext = (...)
% CM.lb_0   = (...)

%% Set meaningful constraints on the design space
% Make a constraint manager object
CM = constraint_manager_cst('Constraint manager 1', SC, SD, SF);    % Create Constraint Manager Object
%  ---- Edit any properties of choice ---- %
    % For example modify the extension factors that quantify how far we can get
    % from the smallest elementary einterval containing the reference designs
    CM.extension_factor_ub_LE = 1.2;                            % Greater value permits thicker airfoils
    CM.extension_factor_ub_M  = 1.2;                            % Greater value permits thicker airfoils
    CM.extension_factor_ub_TE = 1.2;                            % Greater value permits thicker airfoils
        
    CM.extension_factor_lb_LE = 1.2;                            % Greater value permits thinner airfoils
    CM.extension_factor_lb_M  = 1.2;                             % Greater value permits thinner airfoils
    CM.extension_factor_lb_TE = 1.2;                            % Greater value permits thinner airfoils

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
    GM.vectorize = false;                                       % Slower but better for preview. You should always set your case up this way, and once you know it works well, activate vectorization.
    GM.use_reference_population = true;                         % Have a nice start!
    GM.vectorize = false;                                       % Go Multicore!
    GM.options_structure
%  ---- End of property edition ---- %
GM.set_options();                                               % Create option structure according to above defined options
GM.options_structure.PopulationSize=12;
GM.options_structure.Generations=150;
GM.build_problem_structure();                                   % Create problem structure using above options and constraints from CM object
GM.options_structure.MutationFcn= gaoptimset(GM.options_structure, 'MutationFcn', {@mutationadaptfeasible, 0.05});
 
%% Run Optimization
% If you want to run your optimization without the GUI, uncomment the three
% following lines
%   GM.start_genetic_optimizer();
%   results = GM.last_results;
%   save(['results-', datestr(now, 'dd-mm-yyyy')], 'results');

% Fields of interest p
ub0_export = CM.ub0(1:(p_upper.order+p_lower.order));
lb0_export = CM.lb0(1:(p_upper.order+p_lower.order));

ub_ext_export = CM.ub_ext(1:(p_upper.order+p_lower.order));
lb_ext_export = CM.lb_ext(1:(p_upper.order+p_lower.order));


disp('== Original Bounds  ==')
disp( 'Notation : a0_top(LE)  a1_top      a2_top      a3_top      a4_top      a5_top(TE)  a0_bot(LE)  a1_bot      a2_bot       a3_bot      a4_bot     a5_bot(TE)')
disp( '--       : --          --          --          --          --          --          --          --          --           --          --         --        ')
disp(['ub0      : ' num2str(ub0_export)])
disp(['lb0      : ' num2str(lb0_export)])
disp(' ')
disp('== Extended Bounds  ==')
disp( 'Notation : a0_top(LE)  a1_top      a2_top      a3_top      a4_top      a5_top(TE)  a0_bot(LE)  a1_bot      a2_bot       a3_bot      a4_bot     a5_bot(TE)')
disp( '--       : --          --          --          --          --          --          --          --          --           --          --         --        ')
disp(['ub_ext   : ' num2str(ub_ext_export)])
disp(['lb_ext   : ' num2str(lb_ext_export)])
disp(' ')
disp('== Airfoils in Pool ==')
disp( 'Notation                       : a0_top(LE)  a1_top      a2_top      a3_top      a4_top      a5_top(TE)  a0_bot(LE)  a1_bot      a2_bot       a3_bot      a4_bot     a5_bot(TE)')
disp( '--                             : --          --          --          --          --          --          --          --          --           --          --         --        ')
for n_foil = 1:length(CM.names_fitted_bounds)
    str_spaces = make_spaces(30 - length(CM.names_fitted_bounds{n_foil}{1}));
    disp([CM.names_fitted_bounds{n_foil}{1}, str_spaces , ' : ' num2str(CM.ref_pop_array_bounds(n_foil , 1:(p_upper.order+p_lower.order)))])
end

