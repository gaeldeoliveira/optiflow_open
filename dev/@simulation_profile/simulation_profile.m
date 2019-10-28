classdef simulation_profile
    %CLUSTER_PROFILE is a simple class for holding data. For now , we
    % just use it as structure template (non-handle, without methods, but
    % with default values) 
    
    properties
        TA1C 	% 1st Bernstein Polynomial Coefficient for Cf relation
        TA2C    % 2nd Bernstein Polynomial Coefficient for Cf relation
        TA3C    % 3rd Bernstein Polynomial Coefficient for Cf relation
        TA4C    % 4th Bernstein Polynomial Coefficient for Cf relation
        TA5C    % 5th Bernstein Polynomial Coefficient for Cf relation
        TA6C    % 6th Bernstein Polynomial Coefficient for Cf relation
        
        TA1H    % 1st Bernstein Polynomial Coefficient for Hstar relation
        TA2H    % 2nd Bernstein Polynomial Coefficient for Hstar relation
        TA3H    % 3rd Bernstein Polynomial Coefficient for Hstar relation
        TA4H    % 4th Bernstein Polynomial Coefficient for Hstar relation
        TA5H    % 5th Bernstein Polynomial Coefficient for Hstar relation
        TA6H    % 6th Bernstein Polynomial Coefficient for Hstar relation
        
        THMIN   % Lower bound of H intervention region
        THMAX   % Upper bound of H intervention region
        TDCF    % Translation offset of Cf parametrization
    end
    
    methods
        function SPR = simulation_profile()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end