% Clean environment
clear all; close all; clc;

% System Work necessary to start
fs = filesep();                 % Folder separator is OS dependent
addpath([cd fs 'src']);         % Add optimizer source code folder
addpath([cd fs 'user_src']);    % Add optimizer user sources (cost functions, etc) folder
addpath([cd fs 'gui']);         % Add gui source code folder
warning('off','MATLAB:DELETE:FileNotFound');

% Get a database
loaded_datase = load('./data/experimental/EDB18dfdd03-e00b-4717-8af4-6ba3892db2e5.mat');
% Make it present
EDB = loaded_datase.EDB;

% Find out airfoil polar cases
index_airfoil_cases = EDB.find_cases_by_IDkind('airfoil_polar');
% Get them into a cell array
EC_cell_airfoil_cases = EDB.EC_cell(index_airfoil_cases(1:10));

% % Now starts the job handler part

% % Now make the objects that we need to run the cases
% System context
SC = system_context; SC.N_cores = 2; SC.set_context;
% Parametrization objects for upper and lower side
p_upper = parametrization( 'cst_upper'  , 8);
p_lower = parametrization( 'cst_lower'  , 8);
% Shape Definition Object
N_dummy_parameters = 1; SD = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'cst88');

% Now make a general simulation protocol
SP = simulation_protocol('general_base_case' , SC);
SP.target_application='RBINKOLIVEIRA_V2';                      % To use Xfoil as your application set this to 'xfoil'
SP.operation = 'alfa_polar_ref_start';                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
SP.operation_parameters = [-15 20 0.5 0];                      % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15

% SP.Re = 1000000;                                               % Reynolds Number (NaN means inviscid)
% SP1.xtr_top=.99;                                                % Set forced transition on top
% SP1.xtr_bot=.99;                                                % Set forced transition on bottom
% SP1.N_crit=9;                                                   % Set amplification factor of tolmienn schlichting waves
% SP1.save_bl_parameters = 'true';                                % Save and read back boundary layer parameters

% And make a simulation worker
SW = simulation_worker('simul_worker_1', SP, [] , SD,  SC);   % Mariline
SW.app_name = 'RBINKOLIVEIRA_V2';                             % For xfoil write SW_1.app_name = 'xfoil'
SW.parallelized=1;                                            % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
SW.SC.N_cores=2;                                              % If parallelized, over how many cores (not cluster compatible)

% Make a little experiment]
% x = EC_cell_airfoil_cases{1}.airfoil_description.x;
% results_list = analyse_airfoils_list(SW, {x});

% Now make cell lists of airfoil coordinates and polar conditions
x_list                = cell(size(EC_cell_airfoil_cases));
polar_conditions_list = cell(size(EC_cell_airfoil_cases));
for n_airfoil_case = 1:length(EC_cell_airfoil_cases)
    x_list{               n_airfoil_case} = EC_cell_airfoil_cases{n_airfoil_case}.airfoil_description.x;
    polar_conditions_list{n_airfoil_case} = EC_cell_airfoil_cases{n_airfoil_case}.polar_conditions;
end
% results_list = analyse_airfoils_list(SW, x_list);
results_list = SW.analyse_airfoils_list_with_polar_conditions(x_list, polar_conditions_list);

% % Now choose a case, and make comparison for that one (no penalty for non-convergence right not)
sim_results = results_list{1}; exp_case = EC_cell_airfoil_cases{1};

cl_accuracy_index = case_comparator.compare_cl_data(sim_results , exp_case);
cm_accuracy_index = case_comparator.compare_cl_data(sim_results , exp_case);
cd_accuracy_index = case_comparator.compare_cl_data(sim_results , exp_case);















