classdef global_cost_function < handle
    %GLOBAL_COST_FUNCTION Execution description class 
    %   Detailed explanation goes here
    
    properties
        name = 'Cost Function Name';
        description = 'String with Cost Function Description';
        
        % Shape Definition Information
        % shape_definition_handle = {};       % Shape definition handle to 
        
        % Definition of simulations to be executed
        %               Workers are specialized: for each protocol there must be a worker
        % simulation_protocol_list = {};        % List of simulation objects to be executed
        simulation_worker_list = {};          % Function handles to simulation execution module
        
        % Definition of post-processing of executed simulations
        interpretation_function_list = {};
        argument_connectivity_matrix = [1]; %#ok<NBRAK>
        % The argument conectivity matrix describes which, and in which
        % order, experiments are passed to the interpretation functions
        
        % Store last experiment results
        last_experiment_set = {};
        % Store last obtained cost (objective) function array
        last_obj_array      = [];
        % Store last called airfoil coordinates array
        last_x_array        = [];
        
        
        %genotype_storage_handle = {};
        %simul_storage_list = {};
        %cost_function_storage_list = {};
        
        % Vectorization and Paralellism setup
        % vectorized = 0;
        % dedicated_cores = 1;
        
        % Handle to system context object
        % SC
        
    end
    
    methods
        function CF = global_cost_function(name)
            % Constructor method for global_cost_function class objects
            CF.name = name;
        end
        
        function val_obj = interpret_sample(CF , results_sample)
            % Interprets all objectives from experiment results, for a
            % single sample
            
            % Determine number of objectives           
            N_obj  = length(CF.interpretation_function_list);
            val_obj = zeros(1,N_obj);
            
            % Determine number of simulations performed on each sample
            N_sim  = length(CF.simulation_worker_list);
            
            % Cycle for interpretation of each objective
            for n_obj = 1:N_obj
                
                % Build sample_results_list for this objective using the
                % argument_connectivity_matrix which specifies the rules
                % for each objective interpreter.                
                sample_results_list = {};                
                for n_sim = 1:N_sim
                    n_argument = CF.argument_connectivity_matrix(n_sim, n_obj);
                    if n_argument ~= 0
                        sample_results_list{n_argument} = results_sample{n_sim}; %#ok<AGROW>
                    end
                end
                
                % Interpret for Current objective
     %           disp('I am here line 72, gcf')
                val_obj(n_obj) = CF.interpretation_function_list{n_obj}.evaluate(sample_results_list);
            end
            disp(['Cost function value  =  ' , num2str(val_obj)]);
            % Done !
        end
        
        function val_obj_list = interpret_population_results(CF, experiment_matrix_cell)
            % Interprets all experiment results contained in 
            % experiment_matrix_cell cell array, acording to the context of
            % this own CF (global cost function) object
            % Data returned in cell list form.
            
            % Determine how many samples there are
            s = size(experiment_matrix_cell);
            N_samples = s(1);
            
            % Allocate cell list of objective values
            val_obj_list = cell(N_samples, 1);
            
            % Interpret experiment results for each sample
            % % Inserting a parfor here requires that CF be rewritten to
            % avoid to data corruption
            for n_sample = 1:N_samples
                % Extract set of results for this sample
                results_sample = {experiment_matrix_cell{n_sample, :}};
                % Interpret all objective for this sample
                val_obj = CF.interpret_sample(results_sample);
                % Send into list
                val_obj_list{n_sample} = val_obj;
            end
            
        end
        
        function experiment_matrix_cell = run_experiments(CF , x_parameters_list)
            % Preallocate results storage space
            N_sim       = length(CF.simulation_worker_list);
            N_samples   = length(x_parameters_list);
            experiment_matrix_cell = cell(N_samples, N_sim);
            
            % Start by making requested simulations
            for n_sim = 1:N_sim
                % Identify worker for current simulation protocol
                SW = CF.simulation_worker_list{n_sim};
                
                 % Run experiment for set of samples
                 current_sim_results = SW.analyse_airfoils_list(x_parameters_list);
                 
                 % Set data in right format
                 for n_sample = 1:N_samples
                    experiment_matrix_cell{n_sample, n_sim} = current_sim_results{n_sample};
                 end
                    
            end            
            % Finished running experiments
        end
                
        function x_parameters_list = cast_array_into_list(CF , x)
            % Restrcutures the parameter array inputed by the optimizer
            % call into a cell list of samples
            s = size(x);
            N_samples = s(1);
            x_parameters_list = cell(N_samples, 1);
                        
            for n_sample = 1:N_samples
                x_parameters_list{n_sample} = x(n_sample, :);
            end
        end
        
        function val_obj_array = cast_list_into_array(CF, val_obj_list)
            % Restrcutures the list array outputed for example by the
            % method CF.interpret_population_results(...) into a float
            % array adequate for returning to the optimizer            
            val_obj_array = cell2mat(val_obj_list);
           
            % Obsolete code for this purpose:
%             list = val_obj_list
%             N_samples = length(list);
%             N_obj = length(list{1});
%             val_obj_array = zeros(N_samples, N_obj);
%             for n_sample = 1:N_samples
%                 val_obj_array(n_sample,:) = list{n_sample};
%             end
        end
        
        function val_obj_array = evaluate_population(CF , x_parameters_array)
            % Evaluate the cost functions for a whole population. Includes
            % preparing, running and interpreting experiments. This is the
            % most relevant function for use in an optimizer.
            % Accepts vectorized inputs and makes use of parallel computing
            % capabilities.
       
            % Transform array supplied by optimizer callback into cell list of parameters
            x_parameters_list = CF.cast_array_into_list(x_parameters_array);
            % Run all experiments
            experiment_matrix_cell = CF.run_experiments(x_parameters_list);
          %  experiment_matrix_cell
          %  experiment_matrix_cell{1}
            % Interpret all experiments
            val_obj_list = CF.interpret_population_results(experiment_matrix_cell);
            
            % Put experiments into a convenient form for optimizer
            val_obj_array = CF.cast_list_into_array(val_obj_list);
            
            % We are done, just store experiments for eventual
            % postprocessing
            CF.last_experiment_set = experiment_matrix_cell;
            % And cost function values
            CF.last_obj_array      = val_obj_array;
            % And cost function values
            CF.last_x_array        = x_parameters_array;
            
        end
    end
    
end

