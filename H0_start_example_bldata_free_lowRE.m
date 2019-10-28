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
 %clc
 
%% System Work necessary to start
% Start by adding the source folders to the matlab path
fs = filesep();                 % Folder separator is OS dependent
addpath([cd fs 'src']);         % Add optimizer source code folder
addpath([cd fs 'user_src']);    % Add optimizer user sources (cost functions, etc) folder
addpath([cd fs 'gui']);         % Add gui source code folder
warning('off','MATLAB:DELETE:FileNotFound');

% Create System Context Object and Set Context
SC = system_context;
%  ---- Edit any properties of choice, like number of cores ---- %
    SC.N_cores = 4;             % N_cores is actually equivalent to number of Threads. (this was written when cores were single threaded on a Core 2 Duo ULV!)
    
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
    SP1.Re = 500000;                                                % Reynolds Number (NaN means inviscid)
    SP1.xtr_top=.99;                                                % Set forced transition on top
    SP1.xtr_bot=.99;                                                % Set forced transition on bottom
    SP1.N_crit=9;                                                   % Set amplification factor of tolmienn schlichting waves
    SP1.save_bl_parameters = 'true';                                % Save and read back boundary layer parameters
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
CFaero.cost_function_handle = @cost_function_CLCD_at_CL;            % Specify which function to use
CFaero.post_function = @(x) abs(1 / x);                             % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFaero.parameter_list = {'cl_i' , 'w_i'};                           % list of parameters names that control cost function execution
CFaero.parameter_values  = {0.5:0.1:1.2 , [0.05 0.1, 0.1 , 0.1 , 0.1 , 0.15 , 0.2 , 0.2 ]};     % list of parameters values that control cost function execution

% Create a simple structural cost function
CFbldata = cost_function('pressure_coefficient_gradient');          % Initialize it and give it a name
CFbldata.cost_function_handle = @cost_function_RMSdCpdX_at_CL;      % Specify which function to use
CFbldata.post_function = @(x) -abs(1 / x);                          % post_function, applied to result of cost_function_handle, for example to invert values or make negative
CFbldata.parameter_list = {'cl_i' , 'w_i'};                         % list of parameters names that control cost function execution
CFbldata.parameter_values  = {0.5:0.1:1.2 , [0.05 0.1, 0.1 , 0.1 , 0.1 , 0.15 , 0.2 , 0.2 ]};   % list of parameters values that control cost function execution




%% Create the global cost function using the single objectives we defined earlier
% Use the global_cost_function class
CFG = global_cost_function('my_first_one');                     % Initialize it and give it a name
    CFG.simulation_worker_list = {SW_1};                            % Provide list of simulation workers to be executed on each cost function call
    CFG.interpretation_function_list = {CFaero , CFbldata };        % Provide list of cost functions for interpretation of results
    CFG.argument_connectivity_matrix = [1 1];                       % Specify how each cost function uses results from experiments
 
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
    CM.extension_factor_lb_M  = 1.4;                            % Greater value permits thinner airfoils
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
    %GM.vectorize = false;                                       % Slower but better for preview. You should always set your case up this way, and once you know it works well, activate vectorization.
    GM.vectorize = true;                                        % Go Multicore!
    GM.use_reference_population = true;                         % Have a nice start!    
%  ---- End of property edition ---- %
GM.set_options();                                               % Create option structure according to above defined options
GM.options_structure.PopulationSize=240;
GM.options_structure.Generations=80;
GM.build_problem_structure();                                   % Create problem structure using above options and constraints from CM object
GM.options_structure.MutationFcn= gaoptimset(GM.options_structure, 'MutationFcn', {@mutationadaptfeasible, 0.05});

%% Test Cost Function (manual)
disp(['Cost function test call : ' num2str(GM.call_cost_function(GM.CM.ref_pop_array_bounds(1,:)))]);

%% Run Optimization with free transition!
% If you want to run your optimization without the GUI, uncomment the three
% following lines
GM.start_genetic_optimizer();
results = GM.last_results;
save([mfilename(), '-results' , datestr(now, 'yyyymmddTHHMMSS')], 'results');
save([mfilename(), '-fullcase', datestr(now, 'yyyymmddTHHMMSS')]);

% %% Run Optimization again, with forced transition!
% SP1.xtr_top=.10;                                                % Set forced transition on top
% SP1.xtr_bot=.15;                                                % Set forced transition on bottom (far enough to deal with stagnation point bypass!)
% 
% GM.start_genetic_optimizer();
% results_forced = GM.last_results;
% save([mfilename(), '-results-forced-' , datestr(now, 'yyyymmddTHHMMSS')], 'results_forced');
% save([mfilename(), '-fullcase-forced-', datestr(now, 'yyyymmddTHHMMSS')]);
% 
% %% Process a bit
% ld_free         =   1 ./ results_free.fval(  :,1);
% RMSdCpdX_free   = - 1 ./ results_free.fval(  :,2);
% 
% ld_forced       =   1 ./ results_forced.fval(:,1);
% RMSdCpdX_forced = - 1 ./ results_forced.fval(:,2);
% 
% %% Remove a strange point
% ld_free         = ld_free([1,3:end]);
% RMSdCpdX_free   = RMSdCpdX_free([1,3:end]);
% 
% %% Plot Pareto Front
% figure(20)
% plot(results_forced.fval(:,1), results_forced.fval(:,2), 'x'); grid on;
% xlabel('cf_1'); ylabel('cf_2'); 
% 
% 
% figure(20)
% plot(ld_free  , RMSdCpdX_free  , 'x'); hold on;
% plot(ld_forced, RMSdCpdX_forced, 'x'); grid on;
% 
% xlabel('L/D weighted at C_l'); 
% ylabel('L2 norm of dCpdX along chord');
% title('Aerodynamic performance vs Tonal noise forcing');
% legend('Clean leading edge, optimised for free transition' , ...
%        'Triped leading edge, optimised for forced transition');
% 
% %% Loading 
% 
% 
