classdef helper_functions
    %HELPER_FUNCTIONS is a simple class collecting static helper functions
    %for the Dogoro and Kirikou solvers
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                     %
    %   Part of -   A simple specialized 2d Vorticity Equation Solver for %
    %               Actuator Disk Flows (Kirikou-Dogoro Suite)            %
    %                                                                     %
    %   Date    :   June 2014 to March 2017                               %
    %   Author  :   Gael de Oliveira                                      %
    %                                                                     %
    %   License :   Case by case written agreement limited to specific    %
    %               applications. Distribution to any individual or       %
    %               organization requires explicit written agreement from %
    %               original author.                                      %
    %                                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
    end
    
    methods
    end
    
    methods(Static)
        function field_vector = get_field_from_object_cell(object_cell , field_name)
            % Returns an array of field values for the field_name field of
            % the objects collected in the object_cell cell array of
            % objects 
            
            % Allocate returned array 
            field_vector = zeros(size(object_cell));
            % Fill it!
            for  n_object = 1:length(field_vector)
                field_vector(n_object) = object_cell{n_object}.(field_name);
            end 
        end
        
        function case_str = make_case_string(b_in, a2_in, Dx)
            % Returns a case description string, used organized storage of
            % figure .mat data files!
            case_str    = ['case' , '-b_in=' , num2str_fix(b_in, 1,2) , ...
                '-a2_in=' , num2str_fix(a2_in, 2,3) , ...
                '-Dx=', num2str_fix(Dx,2,3)];
        end
        
        function mod_str = num2str_fix(a, n_lead_target, n_trail_target)
            % Convert a NUmber into a String with Fixed Leading and Trailing Zeros
            % a = 1.25;
            % n_lead_target  = 2;
            % n_trail_target = 3;
            
            base_str = num2str(a);
            
            % Take care of leading zeros
            n_lead_base   = find(base_str == '.') - 1; % floor(log10(a)) + 1;
            if isempty(n_lead_base)
                n_lead_base = length(base_str);
            end
            n_lead_to_add = max(n_lead_target - n_lead_base, 0);
            leading_zeros = char(ones(1,n_lead_to_add) * num2str(0));
            
            % Take care of trailing zeros
            n_trail_base = length(base_str) - find(base_str == '.');
            if isempty(n_trail_base)
                n_trail_to_add  = n_trail_target;
                trailing_zeros = char(ones(1,n_trail_to_add) * num2str(0));
                mod_str = [leading_zeros , base_str , '.' , trailing_zeros];
            else
                n_trail_to_add  = n_trail_target - n_trail_base;
                if n_trail_to_add > 0
                    trailing_zeros = char(ones(1,n_trail_to_add) * num2str(0));
                    mod_str = [leading_zeros , base_str , trailing_zeros];
                else
                    mod_str = [leading_zeros , base_str(1:end+n_trail_to_add)];
                end
            end
            % Done! Return!
        end
        
    end
    
end

