classdef mixed_field < handle
    %MIXED_FIELD is a class for storing the mixed momentum field
    %   (u_tilde), generate its mesh, facilitate the computation
    %   derivatives and handling of boundary conditions!
    
    %   Architectural Note:
    %       Separate state from context:
    %           context belongs to objects
    %           state belongs to objects 
    
    properties
        
        %u_tilde         % Mesh of u_tilde values (don't know if this will be maintained, separate state from context, context belongs to objects)
                        % zy coordinates defined in the CM mesh object
        
        CM              % Crossflow Mesh Object Used to Initialize Mixed Field Mesh
        
        
    end
    
    methods
        % Constructor Function
        function MF = mixed_field(CM)
            % Process Inputs
            MF.CM = CM;
        end
        
        function u_tilde = initial_u_tilde(MF)
            % Returns a mesh of u_tilde values initialized to 0, using the
            % CM mesh definition as a reference
            
            % Allocate u_tilde mesh (and set at 0 for initial state)
            u_tilde  = zeros(size(MF.CM.z_mesh));
            % Return mesh as output, don't store to avoid messing state and
            % context (this will help handling different time integration 
            % schemes seamlessly!)
        end
        
        function [u_tilde, du_tilde_dz, du_tilde_dy, lap_u_tilde] = ...
                          filtered_u_tilde_gradients_laplacian(MF, u_tilde)
            % Returns u_tilde, its gradient and laplacian over
            % mesh with filtering for BC enforcement and robustness tricks
            %
            %
            % Returns filtered version of u_tilde including BCs:
            %       -> Periodic BC on sides: u_tilde, gradient and
            %            laplacian set to equal values along primary and 
            %            side type 1 (center of vortex pair) symmetry lines
            %
            % Returns the filtered gradient of u_tilde including BCs:
            %       -> Homogeneous Neumann (fixed gradient) on Wall and Top 
            %       -> Periodic BC on sides: u_tilde, gradient and
            %            laplacian set to equal values along primary and 
            %            side type 1 (center of vortex pair) symmetry lines
            %
            % Returns the filtered laplacian of u_tilde including tricks:
            %       -> No difusion to wall: fudge laplacian to 0 along wall
            %           to avoid transmission of momentum to wall
            %       -> Periodic BC on sides: u_tilde, gradient and
            %            laplacian set to equal values along primary and 
            %            side type 1 (center of vortex pair) symmetry lines
            
            % Get derivatives of u_tilde field
            [du_tilde_dz , du_tilde_dy] = gradient(u_tilde, MF.CM.z_range, MF.CM.y_range);
            % Get Laplacian of u_tilde field
            lap_u_tilde = del2(u_tilde, MF.CM.z_range, MF.CM.y_range);
            
            % Reduce Instability risks at wall (first line, y=0)
            % Tilde Field
            du_tilde_dz(1,:) = 0;
            lap_u_tilde(1,:) = 0;
            
            % Set a boundary condition on top (last line, y=max)
            % Tilde Field
            du_tilde_dz(end,:) = 0;
            du_tilde_dy(end,:) = 0;
            
            lap_u_tilde(end,:) = 0;
            
            % Replicate all fields from center to sides (for periodic BCs)
            % Tilde Field
            j_center = MF.CM.j_center;
            
            u_tilde(:,1  )     = u_tilde(:,j_center);
            u_tilde(:,end)     = u_tilde(:,j_center);
            
            du_tilde_dz(:,1  ) = du_tilde_dz(:,j_center);
            du_tilde_dz(:,end) = du_tilde_dz(:,j_center);
            du_tilde_dy(:,1  ) = du_tilde_dy(:,j_center);
            du_tilde_dy(:,end) = du_tilde_dy(:,j_center);
            
            lap_u_tilde(:,1  ) = lap_u_tilde(:,j_center);
            lap_u_tilde(:,end) = lap_u_tilde(:,j_center);
        end
        
    end
    
end

