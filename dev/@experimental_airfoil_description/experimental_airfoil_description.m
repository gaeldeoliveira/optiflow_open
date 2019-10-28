classdef experimental_airfoil_description
    %AIRFOIL_DESCRIPTION is a simple class for holding data. For
    %now, we just use it as structure template (non-handle, without methods)
    
    properties
        x                   % Shape Description vector parameter (upper cst parameters + lower cst parameters + trailing edge thickness)
        N_dummy_parameters  % Number of dummy parameters (usually 1, for 
        
        tx_coordinates      % Scaled, Translated and Rotated Airfoil Coordinates
        tz_coordinates      % Scaled, Translated and Rotated Airfoil Coordinates
        
        tx_coordinates_raw  % Original Airfoil Coordinates loaded from file
        tz_coordinates_raw  % Original Airfoil Coordinates loaded from file
        
        theta               % Theta rotation (rad) of airfoil, from raw to standard coordinates
        scale_factor        % Factor by which chord is divided, when moving from raw to standard coordinates
        translation_factor  % Offset by which raw airfoil was moved to set LE at origin in standard coordinates
        
        airfoil_filename    % Filename of airfoil file (can include folder name)
        airfoil_name        % Name of airfoil (for display uses only)
    end
    
    methods
        function airfoil_description = experimental_airfoil_description()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end

