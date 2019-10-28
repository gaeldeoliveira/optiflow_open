classdef global_function_wrapper < handle
    %GLOBAL_FUNCTION_WRAPPER Takes care of interfaces with the minimization
    %algorithm for machine learning of closure relations.
    
    properties
        % Inputed fields
        CPD       % Handle to case profile dispatcher object
        SPR_lb    % Describes lower bounds of parameter values (empty for fixed parameters)
        SPR_0     % Describes initial (reference) point (all points filled)
        SPR_ub    % Describes upper bounds of parameter values (empty for fixed parameters)
        SPR_free  % Highlights free (=1, changed directly by optimizer) and fixed (0, paired to something else)
        SPR_index % Describes correspondence to indices in vector passed to optimizer (0 corresponds to copying values from reference case)
        
        % Processed fields
        SPR_fields    % Cell array with list of fields in SPR objects (to allow evolution and future features!)
        N_free_fields % Number of free fiels (that can be changed by minimization algorithm
        
        x_lb      % Lower bounds of parameter vector (used by minimization algorithm)
        x_0       % Initial point parameter vector (used by minimization algorithm)
        x_ub      % Upper bounds of parameter vector (used by minimization algorith
        
        accuracy_profile_0  % Accuracy profile of initial (reference) point
        APE_0               % APE object of initial (reference) point
        % Data-logging fields
        log_requested_scalar_accuracy_metric_array
        log_requested_normalized_accuracy_profile_cell
        log_requested_x_cell
        
    end
    
    methods
        function GFW = global_function_wrapper(CPD, SPR_lb, SPR_0, SPR_ub, SPR_free, SPR_index)
            %GLOBAL_FUNCTION_WRAPPER constructor receives a
            %case_profile_dispatcher object and five simulation_profile
            %objects as input:
            %   CPD       % Handle to case profile dispatcher object
            %   SPR_lb    % Describes lower bounds of parameter values (empty for fixed parameters)
            %   SPR_0     % Describes initial (reference) point (all points filled)
            %   SPR_ub    % Describes upper bounds of parameter values (empty for fixed parameters)
            %   SPR_free  % Highlights free (=1, changed directly by optimizer) and fixed (0, paired to something else)
            %   SPR_index % Describes correspondence to indices in vector passed to optimizer (0 corresponds to copying values from reference case)
            
            % Store Handle to Case Profile Dispatcher
            GFW.CPD         = CPD;
            % Store Simulation Profile Objects
            GFW.SPR_lb      = SPR_lb;
            GFW.SPR_0       = SPR_0;
            GFW.SPR_ub      = SPR_ub;
            GFW.SPR_free    = SPR_free;
            GFW.SPR_index   = SPR_index;
            % Initialize accumulator arrays
            GFW.log_requested_scalar_accuracy_metric_array     = [];
            GFW.log_requested_normalized_accuracy_profile_cell = cell(0);
            GFW.log_requested_x_cell                           = cell(0);
            
            
        end
        
        function process_initial_fields(GFW)
            % Make parameter vector
            GFW.SPR_fields = fieldnames(GFW.SPR_0);
            
            % Find maximum index number
            GFW.N_free_fields = 0;
            for n_field = 1:length(GFW.SPR_fields)
                if GFW.SPR_free.(GFW.SPR_fields{n_field}) == 1
                    GFW.N_free_fields = GFW.N_free_fields + 1;
                end
            end
            
            % Allocate parameter vectors
            GFW.x_0  = zeros(1,GFW.N_free_fields);
            GFW.x_lb = zeros(1,GFW.N_free_fields);
            GFW.x_ub = zeros(1,GFW.N_free_fields);
            
            % And fill them in
            % Loop through all fields
            for n_field = 1:length(GFW.SPR_fields)
                % Only consider the free ones (for making parameter vectors)
                if GFW.SPR_free.(GFW.SPR_fields{n_field}) == 1
                    % Get parameter index of current field
                    i_parameter       = GFW.SPR_index.(GFW.SPR_fields{n_field});
                    % Fill in initial point parameter vector
                    GFW.x_0(i_parameter)  = GFW.SPR_0.(GFW.SPR_fields{n_field});
                    % Fill in lower bounds parameter vector
                    GFW.x_lb(i_parameter) = GFW.SPR_lb.(GFW.SPR_fields{n_field});
                    % Fill in upper bounds parameter vector
                    GFW.x_ub(i_parameter) = GFW.SPR_ub.(GFW.SPR_fields{n_field});
                end
            end
        end
        
        function make_initial_point_accuracy_profile(GFW)
            % Used for computing normalized accuracy_profiles
            [GFW.accuracy_profile_0, GFW.APE_0] = GFW.CPD.dispatch_simulation_profile(GFW.SPR_0);
        end
        
        function make_initial_point_accuracy_profile_on_minibatch(GFW, index_minibatch)
            % Used for computing normalized accuracy_profiles
            % [GFW.accuracy_profile_0, GFW.APE_0] = GFW.CPD.dispatch_simulation_profile(             GFW.SPR_0);
            [GFW.accuracy_profile_0, GFW.APE_0] = GFW.CPD.dispatch_simulation_profile_on_minibatch(GFW.SPR_0, index_minibatch);
        end
        
        function SPR = make_SPR_from_parameters(GFW, x)
            % % And reconstruct an SPR from parameter vector (receive x, return SPR)
            % Create virgin simulation profile
            SPR = simulation_profile();
            % Now fill it in!
            for n_field = 1:length(GFW.SPR_fields)
                % Get parameter index of current field
                i_parameter       = GFW.SPR_index.(GFW.SPR_fields{n_field});
                
                % If field is defined from parameter vector (may be non-free, but duplicated)
                if not(i_parameter == 0)
                    % Fill it in from parameter
                    SPR.(GFW.SPR_fields{n_field}) = x(i_parameter);
                else
                    % Duplicate it from template (initial point)
                    SPR.(GFW.SPR_fields{n_field}) = GFW.SPR_0.(GFW.SPR_fields{n_field});
                end
            end
        end
        
        function log_data(GFW, scalar_accuracy_metric, normalized_accuracy_profile, x)
            % Data logging
            GFW.log_requested_scalar_accuracy_metric_array(end+1)     = scalar_accuracy_metric;
            GFW.log_requested_normalized_accuracy_profile_cell{end+1} = normalized_accuracy_profile;
            GFW.log_requested_x_cell{end+1}                           = x; 
        end
        
        function [scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE] = accuracy_profiles_from_SPR(GFW, SPR)
            % Request the evaluation of a simulation profile
            [accuracy_profile, APE] = GFW.CPD.dispatch_simulation_profile(SPR);
            % Normalize it to reference
            normalized_accuracy_profile = case_comparator.normalize_accuracy_profile(accuracy_profile, GFW.accuracy_profile_0);
            % Make scalar accuracy metric for all three components
            scalar_accuracy_metric = case_comparator.make_scalar_accuracy_metric_from_normalized_accuracy_profile(normalized_accuracy_profile);
            % Make plots
            GFW.make_plots(normalized_accuracy_profile, accuracy_profile)
        end
        
        function [scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE] = accuracy_profiles_from_SPR_on_minibatch(GFW, SPR, index_minibatch)
            % Request the evaluation of a simulation profile
            %[accuracy_profile, APE] = GFW.CPD.dispatch_simulation_profile(SPR);
            [accuracy_profile, APE] = GFW.CPD.dispatch_simulation_profile_on_minibatch(SPR, index_minibatch);
            % Normalize it to reference
            % normalized_accuracy_profile = case_comparator.normalize_accuracy_profile(             accuracy_profile, GFW.accuracy_profile_0);
            normalized_accuracy_profile = case_comparator.normalize_accuracy_profile_on_minibatch(accuracy_profile, GFW.accuracy_profile_0, index_minibatch);
            % Make scalar accuracy metric for all three components
            scalar_accuracy_metric = case_comparator.make_scalar_accuracy_metric_from_normalized_accuracy_profile(normalized_accuracy_profile);
            % Make plots
            GFW.make_plots(normalized_accuracy_profile, accuracy_profile);
        end
        
        function [scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE, SPR] = accuracy_profiles_from_parameters(GFW, x)
            % Reconstruct an SPR from parameter vector (receive x, return SPR)
            SPR = make_SPR_from_parameters(GFW, x);
            % Compute corresponding accuracy profiles
            [scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE] = accuracy_profiles_from_SPR(GFW, SPR); 
            % And log results
            GFW.log_data(scalar_accuracy_metric, normalized_accuracy_profile, x)
        end
        
        function [scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE, SPR] = accuracy_profiles_from_parameters_on_minibatch(GFW, x, index_minibatch)
            % Reconstruct an SPR from parameter vector (receive x, return SPR)
            SPR = make_SPR_from_parameters(GFW, x);
            % Compute corresponding accuracy profiles
            [scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE] = accuracy_profiles_from_SPR_on_minibatch(GFW, SPR, index_minibatch); 
            % And log results
            GFW.log_data(scalar_accuracy_metric, normalized_accuracy_profile, x)
        end
        
        function scalar_accuracy_metric = scalar_metric_from_parameters(GFW, x)
            % Determine number of requested data points !
            N_x                     = size(x,1);
            % Allocate scalar_accuracy_metric output as array !
            scalar_accuracy_metric  = zeros(N_x, 1);
            % Run through requested cases! (and log them as it goes)
            for n_x = 1:N_x
                scalar_accuracy_metric(n_x) = accuracy_profiles_from_parameters(GFW, x(n_x,:));
                lg.msg(GFW,  '--- New simulation profile requested ---');
                lg.msg(GFW, ['                  n_log = ' , num2str(length(GFW.log_requested_x_cell))]);
                lg.msg(GFW, ['                    n_x = ' , num2str(n_x)]);
                lg.msg(GFW, [' scalar_accuracy_metric = ' , num2str(scalar_accuracy_metric(n_x))]);
                lg.msg(GFW, ['                      x = ' , num2str(x(n_x,:))])
            end
        end
        
        function scalar_accuracy_metric = scalar_metric_from_parameters_on_minibatch(GFW, x, index_minibatch)
            % Determine number of requested data points !
            N_x                     = size(x,1);
            % Allocate scalar_accuracy_metric output as array !
            scalar_accuracy_metric  = zeros(N_x, 1);
            % Run through requested cases! (and log them as it goes)
            for n_x = 1:N_x
                scalar_accuracy_metric(n_x) = accuracy_profiles_from_parameters_on_minibatch(GFW, x(n_x,:), index_minibatch);
                lg.msg(GFW,  '--- New simulation profile requested ---');
                lg.msg(GFW, ['                  n_log = ' , num2str(length(GFW.log_requested_x_cell))]);
                lg.msg(GFW, ['                    n_x = ' , num2str(n_x)]);
                lg.msg(GFW, [' scalar_accuracy_metric = ' , num2str(scalar_accuracy_metric(n_x))]);
                lg.msg(GFW, ['                      x = ' , num2str(x(n_x,:))])
                lg.msg(GFW, ['            k_minibatch = ' , num2str(length(index_minibatch))])
                
            end
        end
        
        
    end
        
    methods(Static)
        function make_plots(normalized_accuracy_profile, accuracy_profile)
            figure(101)
            subplot(221)
            hist(accuracy_profile.cl_accuracy_metric_array, 6)
            title(['Global Cl accuracy =' num2str(accuracy_profile.cl_global_accuracy)])
            xlabel('Cl - Lift Coefficient RMS error'); ylabel('Number of Airfoils');
            
            subplot(222)
            hist(accuracy_profile.cm_accuracy_metric_array, 6)
            title(['Global Cm accuracy =' num2str(accuracy_profile.cm_global_accuracy)])
            xlabel('Cm - Pitch Coefficient RMS error'); ylabel('Number of Airfoils');
            
            subplot(223)
            hist(accuracy_profile.cd_accuracy_metric_array, 6)
            title(['Global Cd accuracy =' num2str(accuracy_profile.cd_global_accuracy)])
            xlabel('Cd - Drag Coefficient RMS error'); ylabel('Number of Airfoils');
            
            figure(102)
            subplot(221)
            hist(normalized_accuracy_profile.cl_accuracy_metric_array, 6)
            title(['Global Cl accuracy =' num2str(normalized_accuracy_profile.cl_global_accuracy)])
            xlabel('Cl - Lift Coefficient RMS error'); ylabel('Number of Airfoils');
            
            subplot(222)
            hist(normalized_accuracy_profile.cm_accuracy_metric_array, 6)
            title(['Global Cm accuracy =' num2str(normalized_accuracy_profile.cm_global_accuracy)])
            xlabel('Cm - Pitch Coefficient RMS error'); ylabel('Number of Airfoils');
            
            subplot(223)
            hist(normalized_accuracy_profile.cd_accuracy_metric_array, 6)
            title(['Global Cd accuracy =' num2str(normalized_accuracy_profile.cd_global_accuracy)])
            xlabel('Cd - Drag Coefficient RMS error'); ylabel('Number of Airfoils');
        end    
    end
end

