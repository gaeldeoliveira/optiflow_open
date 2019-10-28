classdef gamultiobj_manager < handle
    %GAMULTIOB_MANAGER Class for managing options and execution of genetic
    % optimizers provided in the Matlab Global Optimization Toolbox
    
    properties
        optimizer_name                                 % gamultiobj or ga
        optimizer_handle
        options_structure
        problem_structure

        fun                                            % Handle to cost function provided to optimizer
        
        CM
        CFG
        ECR
        
        use_eqconstraint_reducer    = false;           % Use the ECR object to reduce number of dimensions of problem using equality constraints
        
        use_ineq_constraints        = false;           % Do not use the inequality constraints by default as matlab seems having trouble to handle them
        use_bounds                  = true;
        use_reference_population    = false;           % If set to true reference population will be inserted into initial population used in GA algorithms
        restricted_initial_range    = false;           % Use restricted upper and lower bounds for creation of initial population
        distance_measure_space      = 'genotype';      % phenotype or genotype
        TolFun                      = 1e-4;
        
        vectorize                   = true;            % Vectorize or not
        
        hybrid_function_postprocessing = false;        % Use fgoalattain to improve pareto front after genetic run
        
        rngstate = [];                                 % If left empty the current state of Random Number Generator will be used, otherwise the value of this field will be used. Setting this is used to replicate resutls
        
        current_population                             % Array containing last requested population (only works for vectorized operation (the object cannot distinguish generations otherwise))
        last_results                                   % Results structure with results from last run
        results_list                                   % Results list storing result structures from all executed runs                
        
        skew_phenotype_space        = false            % Skew phenotype space
        
    end
    
    methods
        function GM = gamultiobj_manager(optimizer_name, CM, CFG, varargin)
            GM.optimizer_name = optimizer_name;
            % Construct optimizer function handle from name
            GM.optimizer_handle = eval([ '@' , GM.optimizer_name]);
            % Generate departing options structure            
            GM.options_structure = gaoptimset(GM.optimizer_name);
            
            % Assign constraint manager handle
            GM.CM = CM;
            % Assig global cost function handle
            GM.CFG = CFG;
            
            if ~isempty(varargin)
                % Assign constraint reducer handle if applicable
                ECR = varargin{1};
                GM.ECR = ECR;
                % Use it as default if it specified
                GM.use_eqconstraint_reducer    = true;
            end
            
            % Create cost function handle
            GM.fun = @(x_array) GM.call_cost_function(x_array);            
        end
        
        function val_obj_array_bis = call_cost_function(GM, x_array)
            
            % Store current set of requested genotypes!
            GM.current_population = x_array;
            
            if GM.use_eqconstraint_reducer == false
                % Call global cost function directly, without any
                % modifications
                val_obj_array = GM.CFG.evaluate_population(x_array);
            end
            
            if GM.use_eqconstraint_reducer == true
                % Here x_array is originally written in the base of the
                % homogeneous solution space of the system Aeq * x = beq
                
                % So, first rewrite into original base
                x_array_original_base = GM.ECR.array_from_null_to_full_space(x_array);
                
                % And now call global cost function
                val_obj_array = GM.CFG.evaluate_population(x_array_original_base);
            end
            
            if GM.skew_phenotype_space == false
                % Do not skew phenotype space (default)
                val_obj_array_bis      =      val_obj_array;
            else                
                % Skew phenotype space
                val_obj_array_bis      = zeros(size(val_obj_array));
                val_obj_array_bis(:,1) = val_obj_array(:,1)     - val_obj_array(:,2) ;
                val_obj_array_bis(:,2) =                          val_obj_array(:,2) ;
            end
            
            % Done !            
        end
        
        function val_obj_array_stripped = call_and_project_cost_function(GM, x_array)
            % Call conventional cost function
            val_obj_array_bis       = GM.call_cost_function(x_array);
            % Strip to first component
            val_obj_array_stripped  = val_obj_array_bis(:,1);
            % Display
            disp([ 'x=' , num2str(x_array(1,:))]);
        end
        
        function set_options(GM)
            %% Set Options
            
            if strcmp(GM.optimizer_name , 'gamultiobj')
                % Plot Pareto Front at end of each generation
                GM.options_structure = gaoptimset(GM.options_structure, 'PlotFcns', @gaplotpareto);
            end
            
            % Output information at each iteration !!! changed by ricardo
            % on 23/02/2014
            %GM.options_structure = gaoptimset(GM.options_structure, 'OutputFcns' , @gaoutputgen);
            GM.options_structure = gaoptimset(GM.options_structure,'Display','iter');
            
            % Force use of linear feasible creation function (insure that
            % lb and ub are respected on first generation)
            GM.options_structure = gaoptimset(GM.options_structure, 'CreationFcn' , @gacreationlinearfeasible);
                    GM.options_structure = gaoptimset(GM.options_structure, 'TolFun' , GM.TolFun);

            % Measure Variety of individuals based on Genotype (instead of default
            % based on phenotype)
            GM.options_structure = gaoptimset(GM.options_structure, 'DistanceMeasureFcn' , {@distancecrowding,GM.distance_measure_space});
            
            GM.options_structure = gaoptimset(GM.options_structure, 'ParetoFraction' , 0.45);
            % Set double mutation rate than default (avoid stalling)
            GM.options_structure = gaoptimset(GM.options_structure, 'MutationFcn', {@mutationadaptfeasible, 0.03}); %%% it used to be 0.02

            if (GM.vectorize == true)
                GM.options_structure = gaoptimset(GM.options_structure, 'Vectorized' , 'on');
            else
                GM.options_structure = gaoptimset(GM.options_structure, 'Vectorized' , 'off');
            end
                
            
            % Insert initial population elements if applicable
            if GM.use_reference_population == true
                GM.options_structure = ...
                    gaoptimset(GM.options_structure, 'InitialPopulation' , GM.CM.ref_pop_array_bounds);
            end
            
            % Set initial range for population
            % Build initial range matrix according to mode selected in
            % options of this object
            if GM.restricted_initial_range == true
                PopInitRange = [GM.CM.lb0 ; GM.CM.ub0];
                GM.options_structure = gaoptimset(GM.options_structure, 'PopInitRange' , PopInitRange);
%             else
%                 PopInitRange = [min(GM.CM.lb_ext) ; max(GM.CM.ub_ext)];
            end
            
            % If applicable, run hybrid function optimization
            if (GM.hybrid_function_postprocessing == true)
               GM.options_structure = gaoptimset(GM.options_structure ,'HybridFcn',@fgoalattain);               
            end

                                    
            
        end
        
        function build_problem_structure(GM)
            % This function creates the problem structure for the current
            % single or multiobjective optimization problem

            problem.fitnessfcn = GM.fun;
            
            % For the original problem, without dimension reduction with 
            % the equality constraints            
            if GM.use_eqconstraint_reducer == false
                problem.nvars = GM.CM.nvars;
                if (GM.use_ineq_constraints == true)
                    problem.Aineq = GM.CM.A_ineq;
                    problem.bineq = GM.CM.b_ineq;
                else
                    problem.Aineq = [];
                    problem.bineq = [];
                end
                problem.Aeq = GM.CM.A_eq;
                problem.beq = GM.CM.b_eq;
                if (GM.use_bounds == true)
                    problem.lb = GM.CM.lb_ext;
                    problem.ub = GM.CM.ub_ext;
                else
                    problem.lb = [];
                    problem.ub = [];
                end
                
                % And set adequate starting population
                if GM.use_reference_population == true
                    GM.options_structure.InitialPopulation = GM.CM.ref_pop_array_bounds;
                end
                
            end
            
            % For the reduced problem, with less dimensions thanks to the
            % use of the equality constraints            
            if GM.use_eqconstraint_reducer == true
                problem.nvars = GM.ECR.nvars_nspace;
                
                % Set all linear constraints to zero, as the inequality
                % ones cannot be considered, and the equality constraints
                % are always satisfied here.
                problem.Aineq = [];
                problem.bineq = [];
                problem.Aeq = [];
                problem.beq = [];                
                
                % Now set bounds!
                problem.lb = GM.ECR.lb_ext_nspace;
                problem.ub = GM.ECR.ub_ext_nspace;
                
                % Add suction constraints if applicable
                if GM.CM.free_suction_area_constraints == true
                    problem.lb(end-1:end) = GM.CM.lb_ext(end-1:end);
                    problem.ub(end-1:end) = GM.CM.ub_ext(end-1:end);
                end
                
                % And set adequate starting population
                if GM.use_reference_population == true
                    GM.options_structure.InitialPopulation = GM.ECR.x_list_fitted;
                end
                
                % And we are done!
                
            end
            
            problem.options = GM.options_structure;
            problem.solver = GM.optimizer_name;
           
           % If we are using the single objective ga, set non linear
           % constraints field (to null as we are not using it currently)
           if strcmp(GM.optimizer_name, 'ga')
               problem.nonlcon = [];
           end

          % Now assign state to random number generator (use specified
          % state if available or take current if user did not specify anything)
%           if ~isempty(GM.rngstate)            
%             problem.rngstate = GM.rngstate;   
%           end
          
          
          % Done! Just store problem structure!
          GM.problem_structure = problem;
          
%           A_ineq = CM.A_ineq
% b_ineq = CM.b_ineq
%           
        end
        
        function start_genetic_optimizer(GM)
            % Run Optimization
            %disp(' lalal alla alalal alalla lala al')
            %GM
            [x,fval,exitflag,output,population,scores] = GM.optimizer_handle(GM.problem_structure);
            
            % Store results!
            results.x = x;
            results.fval = fval;
            results.exitflag = exitflag;
            results.output = output;
            results.population = population;
            results.scores = scores;
%            results.rng_start_state = GM.problem.rngstate;
            
            % Set these results as the last ones
            GM.last_results = results;
            % Add results to result series list (for multiple runs)
            GM.results_list{length(GM.results_list)+1} = GM.last_results;
        end
        
        function rerun_genetic_optimizer(GM)
%             % Restart with results from last optimization                
%             GM.options_structure.InitialPopulation = GM.last_results.population;
%             % Update problem structure with new options
%             GM.build_problem_structure();
            
            GM.problem_structure.options.InitialPopulation = GM.last_results.population;
            
            % Run Optimization
            [x,fval,exitflag,output,population,scores] = GM.optimizer_handle(GM.problem_structure);
            
            % Store results!
            results.x = x;
            results.fval = fval;
            results.exitflag = exitflag;
            results.output = output;
            results.population = population;
            results.scores = scores;
 %           results.rng_start_state = GM.problem.rngstate;
            
            % Set these results as the last ones
            GM.last_results = results;
            % Add results to result series list (for multiple runs)
            GM.results_list{length(GM.results_list)+1} = GM.last_results;                        
        end
            
    end
    
end

