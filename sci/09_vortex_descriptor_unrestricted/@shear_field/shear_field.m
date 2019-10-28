classdef shear_field < handle
    %SHEAR_FIELD is a class for generating the pure shear momentum field
    %   (u_bar) over a crossflow mesh (CM), computing its derivatives and
    %   applying various fudging tricks for improved robustness!
    
    properties
        
        CM              % Crossflow Mesh Object Used to Initialize Mixed Field Mesh
        
        SP              % Swafford Profile Object
        
    end
    
    methods
        % Constructor Function
        function SF = shear_field(CM)
            % Store Inputs
            SF.CM = CM;
            % Create Swafford Profile Object
            SF.SP = swafford_profile();
            
        end
        
        function u_bar_at_y = u_bar_at_y(SF, hk, rt, ue, nue, y)
            % Computes swafford profile at height y, which can be a scalar,
            % 1d-array or 2d-mesh. Scalar arguments expected:
            %       hk          - (kinematic) shape factor
            %       rt          - Reynolds theta (momentum thickness Reynolds number)
            %       ue          - Edge velocity (dimensional or not!)
            %       ue_over_nue - Ratio of edge velocity to kinematic
            %                           viscosity  (for Reynolds number!)
            %       y           - Distance to wall (dimensional, [m])
            %       
            
            % Compute ue_over_nue
            ue_over_nue = ue / nue;
            
            % Generate Scaled Y coordine over mesh (call static method)
            y_over_theta = shear_field.y_over_theta(rt, ue_over_nue, y);
            
            % Compute u_over_ue from Swafford profile
            u_over_ue    = SF.u_over_ue_scaled_y(hk, rt, y_over_theta);
            
            % Compute u_bar from u_over_ue_mesh, that is rescale into
            % suitable velocity coordinates
            u_bar_at_y   = u_over_ue * ue;
            
        end
        
        function u_bar = u_bar_over_mesh(SF, hk, rt, ue, ue_over_nue)
            % Computes swafford profile over mesh define in crossflow_mesh
            % (CM) object. Scalar arguments expected:
            %       hk          - (kinematic) shape factor
            %       rt          - Reynolds theta (momentum thickness Reynolds number)
            %       ue          - Edge velocity (dimensional or not!)
            %       ue_over_nue - Ratio of edge velocity to kinematic
            %                           viscosity  (for Reynolds number!)
            
            % Generate Scaled Y coordine over mesh (call static method)
            y_over_theta_mesh = shear_field.y_over_theta(rt, ue_over_nue, SF.CM.y_mesh);
            
            % Compute u_over_ue from Swafford profile
            u_over_ue_mesh = SF.u_over_ue_scaled_y(hk, rt, y_over_theta_mesh);
            
            % Compute u_bar from u_over_ue_mesh, that is rescale into
            % suitable velocity coordinates
            u_bar          = u_over_ue_mesh * ue;
        end
          
        function u_over_ue = u_over_ue_scaled_y(SF, hk, rt, y_over_theta)
            % Wrapper for swafford profile evaluator: hk and rt must be
            % scalars while y_over_theta can be an array (nD I believe,
            % only tested with scalar, 1d array and 2d mesh)
            u_over_ue = SF.SP.u_over_ue(hk, rt, y_over_theta);
        end
        
        
        function [u_bar, du_bar_dz, du_bar_dy, lap_u_bar] = ...
                                filtered_u_bar_gradients_laplacian(SF, u_bar)
            
            
            % Get Gradient of Bar Field
            [du_bar_dz , du_bar_dy] = gradient(u_bar, SF.CM.z_range, SF.CM.y_range);
            % Get Laplacian of Bar Field
            lap_u_bar = del2(u_bar, SF.CM.z_range, SF.CM.y_range);
            
            % Reduce Instability risks at wall (first line, y=0)
            % Bar Field
            du_bar_dz(1,:) = 0;
            lap_u_bar(1,:) = 0;
            
            % Set a boundary condition on top
            % Bar Field
            du_bar_dz(end,:) = 0;
            du_bar_dy(end,:) = 0;
            
            lap_u_bar(end,:) = 0;
            
            % Replicate fields from center to sides (for periodic BCs)
            % Bar Field
            j_center = SF.CM.j_center;
            
            u_bar(:,1  )       = u_bar(:,j_center);
            u_bar(:,end)       = u_bar(:,j_center);
            
            du_bar_dz(:,1  )   = du_bar_dz(:,j_center);
            du_bar_dz(:,end)   = du_bar_dz(:,j_center);
            du_bar_dy(:,1  )   = du_bar_dy(:,j_center);
            du_bar_dy(:,end)   = du_bar_dy(:,j_center);
            
            lap_u_bar(:,1  )   = lap_u_bar(:,j_center);
            lap_u_bar(:,end)   = lap_u_bar(:,j_center);
          
        end
        
    end
    
    methods(Static)
        function y_over_theta = y_over_theta(rt, ue_over_nue, y)
            % Computes y_over_theta from y, rt (Reynolds Theta) and
            % ue_over_nue. Works for scalar rt, scalar ue_over_nue and nD
            % arrays of y
            
            % Compute momentum thickness
            theta = rt / ue_over_nue;
            % Divide y coordinate by momentum thickness
            y_over_theta = y ./ theta;
            
        end
        
    end
    
end

