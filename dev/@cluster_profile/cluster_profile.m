classdef cluster_profile
    %CLUSTER_PROFILE is a simple class for holding data. For now , we
    % just use it as structure template (non-handle, without methods, but
    % with default values) 
    
    properties
        n_cores_per_node = 2        % Number of cores per node
    end
    
    methods
        function CPR = cluster_profile()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end