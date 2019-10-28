classdef airfoil_polar_evaluator < handle
    %AIRFOIL_POLAR_EVALUATOR is a handle class meant to support the
    % evaluation of a collection of airfoil polars and its comparison with
    % experimental results. This is a handle class meant to be instanciated
    % on each local node (not core, the evaluator uses simulation work
    % facilities to run polars in parallel over the cores of the node).
    
    properties
        CPR          % Cluser Profile determining the type of parallelism expected from the APE object
        SPR          % Simulation Profile used to taylor the behaviour the simulation code
        EC_cell      % A cell array of experimental cases with IDkind 'airfoil_polar' (no runtime check)
        u_order = 8  % Order of upper parametrization (can be edited before initialize_local_context())
        l_order = 8  % Order of lower parametrization (can be edited before initialize_local_context())
        
        SC           % (local) System context object
        p_upper      % (local) Upper side cst parametrization object
        p_lower      % (local) Lower side cst parametrization object
        SD           % (local) Shape definition object
        SP           % (local) Simulation protocol object
        SW           % (local) Simulation worker object

        results_list % (value) Cell array of simulation results
        cl_accuracy_metric_array % Array of Cl accuracy metrics of each simulated airfoil
        cm_accuracy_metric_array % Array of Cm accuracy metrics of each simulated airfoil
        cd_accuracy_metric_array % Array of Cd accuracy metrics of each simulated airfoil
    end
    
    methods
        function APE = airfoil_polar_evaluator(CPR, SPR, EC_cell)
            %AIRFOIL_POLAR_EVALUATOR Constructor method receives:
            %   EC_cell -  A cell array of experimental cases with IDkind
            %              'airfoil_polar' (no runtime check)
            %   SPR     -  A simulation profile used to taylor the
            %               behaviour the simulation code
            %   CPR     -  A cluster profile determining the type of
            %               parallelism expected from the APE object
            
            % Store received values
            APE.CPR     = CPR    ;
            APE.SPR     = SPR    ;
            APE.EC_cell = EC_cell;
        end
        
        function initialize_local_context(APE)
            % % Now make the objects that we need to run the cases
            % System context
            APE.SC = system_context; 
            APE.SC.N_cores = APE.CPR.n_cores_per_node; 
            APE.SC.set_context;
            % Parametrization objects for upper and lower side
            APE.p_upper = parametrization( 'cst_upper'  , APE.u_order);
            APE.p_lower = parametrization( 'cst_lower'  , APE.l_order);
            % Shape Definition Object
            N_dummy_parameters = 1; 
            APE.SD = shape_definition_cst(APE.p_upper , APE.p_lower , N_dummy_parameters, 'cst88');
            
            % Now make a general simulation protocol
            APE.SP = simulation_protocol('ape_ml_case' , APE.SC);
            APE.SP.target_application   = 'RBINKOLIVEIRA_V2';                  % To use Xfoil as your application set this to 'xfoil'
            APE.SP.operation            = 'alfa_polar_ref_start';              % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice
            APE.SP.operation_parameters = [-15 20 0.25 0];                     % Make polar from -15 to 20 degrees in 0.25 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
            APE.SP.SPR                  = APE.SPR;                             % Write 
            APE.SP.EPSV                 = 1e-6;                                % Set more stringent convergence requirements
            
            % And make a simulation worker
            APE.SW = simulation_worker('simul_worker_1', APE.SP, [] , APE.SD,  APE.SC);    % Mariline
            APE.SW.app_name = 'RBINKOLIVEIRA_V2';                              % For xfoil write SW_1.app_name = 'xfoil'
            APE.SW.parallelized=1;                                             % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
            APE.SW.SC.N_cores=APE.CPR.n_cores_per_node;                        % If parallelized, over how many cores (not cluster compatible)
        end
        
        function run_simulations(APE)                
            % Allocate cell lists of airfoil coordinates and polar conditions
            x_list                = cell(size(APE.EC_cell));
            polar_conditions_list = cell(size(APE.EC_cell));
            % Fill them in!
            for n_airfoil_case = 1:length(APE.EC_cell)
                x_list{               n_airfoil_case} = APE.EC_cell{n_airfoil_case}.airfoil_description.x;
                polar_conditions_list{n_airfoil_case} = APE.EC_cell{n_airfoil_case}.polar_conditions;
            end
            % Run simulations!
            APE.results_list = APE.SW.analyse_airfoils_list_with_polar_conditions(x_list, polar_conditions_list);
        end
        
        function [cl_global_accuracy , cm_global_accuracy , cd_global_accuracy] = compare_simulation_results(APE)
            % Allocate arrays for accuracy indexes
            APE.cl_accuracy_metric_array = zeros(size(APE.EC_cell));
            APE.cm_accuracy_metric_array = zeros(size(APE.EC_cell));
            APE.cd_accuracy_metric_array = zeros(size(APE.EC_cell));
            % Compare each result!
            for n_EC = 1:length(APE.EC_cell)
                try
                sim_results = APE.results_list{n_EC}; exp_case = APE.EC_cell{n_EC};
                APE.cl_accuracy_metric_array(n_EC) = case_comparator.compare_cl_data(sim_results , exp_case);
                APE.cm_accuracy_metric_array(n_EC) = case_comparator.compare_cm_data(sim_results , exp_case);
                APE.cd_accuracy_metric_array(n_EC) = case_comparator.compare_cd_data(sim_results , exp_case);
                catch
                    disp(n_EC);
                    disp('Screwed up');
                   
                end
            end
            % Group the cases!
            cl_global_accuracy = case_comparator.group_case(APE.cl_accuracy_metric_array);
            cm_global_accuracy = case_comparator.group_case(APE.cm_accuracy_metric_array);
            cd_global_accuracy = case_comparator.group_case(APE.cd_accuracy_metric_array);
        end
        
        function [ap , experiment_result]= make_polar_on_airfoil(APE, x_airfoil, polar_conditions)
            % Run experiment
            experiment_results       = APE.SW.analyse_airfoils_list_with_polar_conditions({x_airfoil}, {polar_conditions});
            % Extract experiment results
            experiment_result        = experiment_results{1};
            % Extract polar
            ap = experiment_result.ap;
        end
    end
end

