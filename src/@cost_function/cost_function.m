classdef cost_function < handle
    %COST_FUNCTION is a handle class containing the description of 
    %   Detailed explanation goes here
    
    properties
        name = {'risoea1'}
        description = {'Risoe A series style cl/cd function'}
        cost_function_handle = @cost_function_risoea_CLCD;
        post_function = @(x) 1 / x;                             % post_function, applied to result of cost_function_handle, for example to invert values or make negative

        parameter_list = {'alpha_i' , 'w_i'};
        parameter_values  = {1:5 , [0.2 0.2 0.2 0.2 0.2]};
        parameter_structure
        
        failure_value = +Inf;       % Value at which to set result in case operations yielded non-finite(NaN, +Inf or -Inf) value
        
        % comand_list = {}
        % import_from_file = 
    end
    
    methods
        function CF = cost_function(name)
            CF.name = name;
        end
        
        function make_parameter_structure(CF)
        % This function builds a parameter structure passed as first
        % argument to the cost function called by cost_function_handle. 
        % The parameters structure fields and values are contained in the
        % parameter_list and parameter_values cell arrays.
            if isempty(CF.parameter_list)
                CF.parameter_structure = [];
            else
                N_fields = length(CF.parameter_list);
                for n_field = 1:N_fields
                    field_name  = CF.parameter_list{n_field};
                    field_value = CF.parameter_values{n_field};
                    CF.parameter_structure.(field_name) = field_value;
                end
            end
        end
        
        function val = evaluate(CF, experiments_results)
            % Evaluates cost function!
                        
            % Remake Parameter Structure on each call to ensure it is up to
            % date
            CF.make_parameter_structure();
            % Call cost function specified in handle (parameters first,
            % experiment results second)
            cval = CF.cost_function_handle(CF.parameter_structure, experiments_results);
            
            % Apply post_function, for example to invert values or make
            % negative
            val = CF.post_function(cval);
            
            % Insure a value compatible with optimizer standards is
            % outputed in case of failure
            if ~isfinite(val)
                val = CF.failure_value;
            end
%            disp(['Cost function value  =  ' , num2str(val)])
        end
    end
    
end

