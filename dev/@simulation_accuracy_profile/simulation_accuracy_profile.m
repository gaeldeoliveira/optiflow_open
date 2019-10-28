classdef simulation_accuracy_profile
    %ACCURACY_PROFILE is a simple class for holding data. It describes the
    % accuracy of a simulation profile. For now, we just use it as structure
    % template (non-handle, without methods). Each SPR is associated with an
    % accuracy profile
    
    
    properties
        cl_global_accuracy          % Accuracy of Cl computation on airfoil_polar cases
        cm_global_accuracy          % Accuracy of Cm computation on airfoil_polar cases
        cd_global_accuracy          % Accuracy of Cd computation on airfoil_polar cases
        
        cl_accuracy_metric_array    % Accuracy of Cl computation for each airfoil_polar case
        cm_accuracy_metric_array    % Accuracy of Cm computation for each airfoil_polar case
        cd_accuracy_metric_array    % Accuracy of Cd computation for each airfoil_polar case
    end
    
    methods
        function accuracy_profile = simulation_accuracy_profile()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end

