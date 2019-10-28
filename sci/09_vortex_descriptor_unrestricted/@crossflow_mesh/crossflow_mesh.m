classdef crossflow_mesh < handle
    %CROSSFLOW_MESH is a handle class that stores a structured mesh of y-z
    % points for discretizing mixing (u_tilde) and shear (u_bar) fields
    % over the crossflow plane.
    %   z - corresponds to the spanwise direction
    %   y - corresponds to the wall normal direction
    % Fielnames are self-explanatory 
    
    
    properties
        % Inputs: Mesh Properties
        n_range_z = 101;%101;    % Number of Elements in z direction (must be odd)
        n_range_y = 101;%101;    % Number of Elements in z direction (odd or not!)
        
        S                        % Half Spacing between vortex pair centers (z direction)
        delta                    % Boundary Layer Height Estimation (y direction)
        
        N_delta = 4;             % Mesh height in delta units (traditional 4)
        Exp_y   = 1.3;           % Lumping exponent for y mesh 
        
        % Outputs: 1-D Point Arrays (Row)
        z_range                  % Range of steps in the spanwise (z) direction
        y_range                  % Range of steps in the normal   (y) direction
        
        % Outputs: 2-D Point Arrays (Mesh)
        z_mesh                   % Mesh of y positions for velocity fields
        y_mesh                   % Mesh of z positions for velocity fields
        
        % Outputs: Symmetry Indices
        j_center                 % Single integer index for z=0 position
        
        % Outputs: Window Width (for integration)
        z_lenght
        
    end
    
    methods
        % Constructor Method
        function CM = crossflow_mesh(S, delta)
            % Store Inputs
            CM.S     = S;
            CM.delta = delta;
            
            % Update Dependent Fields
            CM.update_dependent_fields();
        end
        
        function make_ranges(CM)
            % Define Steps in Spanwise Direction (constant spacing)
            CM.z_range =   CM.S     *  linspace(-2, 2, CM.n_range_z)        ;
            % Define Steps in Normal Direction (constant spacing)
            CM.y_range = CM.N_delta*CM.delta * (linspace( 0, 1, CM.n_range_y)).^(CM.Exp_y);  % 4*CM.delta * (linspace( 0, 1, CM.n_range_y)).^(1.3)
            % Make lenght of z window, for integration of thicknesses
            CM.z_lenght = max(CM.z_range) - min(CM.z_range);
        end
        
        function make_meshes(CM)
            % Make meshes
            [CM.z_mesh, CM.y_mesh] = meshgrid(CM.z_range, CM.y_range);
        end
        
        function make_symmetry_indices(CM)
            % Make Symmetry Indices 
            
            % Integer index of points with z=0 position (in z_range array
            % and second index of mesh arrays!)
            CM.j_center = find(CM.z_range == 0); 
            
            
        end
        
        function update_dependent_fields(CM)
            % Update all fields upon insertion of new values!
            
            % Make ranges
            CM.make_ranges();
            % Make Meshes
            CM.make_meshes();
            % Make Symmetry Indices
            CM.make_symmetry_indices();
        end
        
    end
    
end

