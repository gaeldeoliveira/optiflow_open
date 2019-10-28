classdef experimental_processed_polars
    %EXPERIMENTAL_PROCESSED_POLARS is a simple class for holding data. For
    %now, we just use it as structure template (non-handle, without methods)
    
    
    properties
        % Returns a structure, processed_polars, with fields:
        % First : alpha, cl polar
        al_alcl_polar    % Angle of attack (deg)
        cl_alcl_polar    % Lift coefficient
        % Second: cl, cd polar
        cl_clcd_polar    % Lift coefficient
        cd_clcd_polar    % Drag coefficient
        % Third : alpha, cl, cd polar
        al_alclcd_polar  % Angle of attack (deg)
        cl_alclcd_polar  % Lift coefficient
        cd_alclcd_polar  % Drag coefficent
        % Fourth: alpha, cm polar
        al_alcm_polar    % Angle of attack (deg)
        cm_alcm_polar    % Picth coefficient
        % Metadata:
        datataset_filnename
        source
        
        
    end
    
    methods
        function EPP = experimental_processed_polars()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end

