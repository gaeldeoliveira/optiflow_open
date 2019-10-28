classdef experimental_bl_run
    %EXPERIMENTAL_PROCESSED_POLARS is a simple class for holding data. For
    %now, we just use it as structure template (non-handle, without methods)
    
    
    properties
        % Stores protected information on 
        x               % [m     ] - Array of x positions at which B.L. was measured
        Ue              % [m/s   ] - Array of edge velocity positions at measurement points 
        theta           % [m     ] - Array of BL momentum displacements at measurement points 
        H               % [adim. ] - Array of shape factors at measurement points 
        
        U_ref           % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
        U_ref_over_nu   % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
        source          % [string] - Identifier for data origin
        
        % Metadata:
        datataset_filnename
        
        
        
    end
    
    methods
        function bl_run = experimental_bl_run()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end