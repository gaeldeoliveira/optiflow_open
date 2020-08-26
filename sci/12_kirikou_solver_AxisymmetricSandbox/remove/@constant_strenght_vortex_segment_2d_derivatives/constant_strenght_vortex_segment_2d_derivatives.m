classdef constant_strenght_vortex_segment_2d_derivatives < handle
    %constant_strenght_vortex_segment
    %   is a class describing a single vortex segment with constant
    %   strenght
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
        % % Independent fields
        x_start             % Vortex segment start coordinates
        y_start             % 
        
        x_end               % Vortex segment end coodinates
        y_end               %
        
        gamma               % Vortex segment strenght
        
        % % Dependent fields
        n_segments          % Number of segments
        
        x_center            % Vortex segment center coordinates
        y_center            %
        
        l                   % Vortex segment lenght
        
        x_l_unit_vector     % Directed Axis Unit vector components
        y_l_unit_vector     %
        
        x_n_unit_vector     % Normal Unit vector components
        y_n_unit_vector     % 
        
        % % Derivatives of Dependent fields
        x_center_dx_end     % Derivatives of Center point coordinates to edges
        x_center_dx_start
        y_center_dy_end
        y_center_dy_start
        
        l_dx_end            % Derivatives of lenght to edge coordinates
        l_dx_start
        l_dy_end
        l_dy_start
        
                            % Derivatives of tangential unit vector to edge coordinates
        x_l_unit_vector_dx_end  
        x_l_unit_vector_dx_start
        y_l_unit_vector_dy_end
        y_l_unit_vector_dy_start
        
                            % Derivatives of normal unit vector to edge coordinates
        x_n_unit_vector_dy_end
        x_n_unit_vector_dy_start
        y_n_unit_vector_dx_end
        y_n_unit_vector_dx_start
    end
    
    methods
        function VS = constant_strenght_vortex_segment_2d_derivatives(x_start, y_start, x_end, y_end, gamma)
            % Constructor Function
            
            % Use classical synchronization routine
            update_independent_fields(VS, x_start, y_start, x_end, y_end, gamma)
        end
        
        function update_independent_fields(VS, x_start, y_start, x_end, y_end, gamma)
            % Updates the independent fields and synchornizes the dependent
            % fields afterwards            
            VS.x_start = x_start(:);
            VS.y_start = y_start(:);
            VS.x_end   = x_end(:);
            VS.y_end   = y_end(:);
            VS.gamma   = gamma(:);
            % Synchronize
            synchronize_dependent_fields(VS)
        end
        
        function synchronize_dependent_fields(VS)
            % Make the dependent fields consistent with the independent
            % fields
            
            VS.n_segments = length(VS.x_end);
            
            VS.x_center = 0.5 * (VS.x_end + VS.x_start);
            VS.y_center = 0.5 * (VS.y_end + VS.y_start);
            
            VS.l = sqrt((VS.x_end-VS.x_start).^2 + (VS.y_end-VS.y_start).^2);
            
            VS.x_l_unit_vector = (VS.x_end - VS.x_start) ./ VS.l;
            VS.y_l_unit_vector = (VS.y_end - VS.y_start) ./ VS.l;
            
            VS.x_n_unit_vector = - VS.y_l_unit_vector;
            VS.y_n_unit_vector =   VS.x_l_unit_vector;
            
            % Synchronize derivatives of the dependent fields
            synchronize_dependent_fields_derivatives(VS)
        end
        
        function synchronize_dependent_fields_derivatives(VS)
            % Make the dependent fields derivatives consistent with the 
            % dependent/independent fields
            
            % Derivatives of Center point coordinates to edges
            VS.x_center_dx_end   = 0.5 * (1                    );
            VS.x_center_dx_start = 0.5 * (         + 1         );
            
            VS.y_center_dy_end   = 0.5 * (1                    );
            VS.y_center_dy_start = 0.5 * (         + 1         );
            
            
            % Derivatives of lenght to edge coordinates
            VS.l_dx_end   = (0.5 ./ VS.l) .*  (  2*(VS.x_end-VS.x_start));
            VS.l_dx_start = (0.5 ./ VS.l) .*  ( -2*(VS.x_end-VS.x_start));
            VS.l_dy_end   = (0.5 ./ VS.l) .*  (  2*(VS.y_end-VS.y_start));
            VS.l_dy_start = (0.5 ./ VS.l) .*  ( -2*(VS.y_end-VS.y_start));
            
            
            % Derivatives of tangential unit vector to edge coordinates
            VS.x_l_unit_vector_dx_end   = - (VS.x_end - VS.x_start) ./ ((VS.l).^2) .* VS.l_dx_end   + (1    ) ./ VS.l;
            VS.x_l_unit_vector_dx_start = - (VS.x_end - VS.x_start) ./ ((VS.l).^2) .* VS.l_dx_start + (  - 1) ./ VS.l;
            
            VS.y_l_unit_vector_dy_end   = - (VS.y_end - VS.y_start) ./ ((VS.l).^2) .* VS.l_dy_end   + (1    ) ./ VS.l;
            VS.y_l_unit_vector_dy_start = - (VS.y_end - VS.y_start) ./ ((VS.l).^2) .* VS.l_dy_start + (  - 1) ./ VS.l;
            
            
            % Derivatives of normal unit vector to edge coordinates
            VS.x_n_unit_vector_dy_end   = - VS.y_l_unit_vector_dy_end;
            VS.x_n_unit_vector_dy_start = - VS.y_l_unit_vector_dy_start;
            
            VS.y_n_unit_vector_dx_end   =   VS.x_l_unit_vector_dx_end;
            VS.y_n_unit_vector_dx_start =   VS.x_l_unit_vector_dx_start;
        end
        
        function [pot , u , v, u_dx, u_dy, v_dx, v_dy] = influence_on_single_point(VS, x, y)
            % Compute induced speed on single point
            [pot , u , v, u_dx, u_dy, v_dx, v_dy] = ...
                constant_strenght_vortex_segment_element_derivatives( ...
                                                         x         , y, ...
                                                         VS.x_start, VS.y_start, ...
                                                         VS.x_end  , VS.y_end, ...
                                                         VS.gamma); 
        end
        
        function [pot , u , v, u_dx, u_dy, v_dx, v_dy] = induced_speed_on_single_point(VS, x, y)
            % Compute induced speed on single point
            
            % Get Single Point Influence
            [pot , u , v, u_dx, u_dy, v_dx, v_dy] = VS.influence_on_single_point(x, y);
            
            % Sum it up
            pot  = sum(pot);
            u    = sum(u);
            v    = sum(v);
            u_dx = sum(u_dx);
            u_dy = sum(u_dy);
            v_dx = sum(v_dx);
            v_dy = sum(v_dy);
        end
        
        function [pot , u , v, u_dx, u_dy, v_dx, v_dy] = influence_on_many_points(VS, x, y)
            disp(['WARNING: VS::influence_on_many_points() was called but it is dorment!'])
%             n_points = size(x);
% 
%             pot = zeros(VS.n_segments, n_points(1), n_points(2));
%             u   = zeros(VS.n_segments, n_points(1), n_points(2));
%             v   = zeros(VS.n_segments, n_points(1), n_points(2));
%           
%             for n_point1 = 1:n_points(1)
%                 for n_point2 = 1:n_points(2)
%                     [pot(:,n_point1 , n_point2) , u(:,n_point1 , n_point2) , v(:,n_point1 , n_point2)] = ...
%                         influence_on_single_point(VS, x(n_point1 , n_point2), y(n_point1 , n_point2));
%                 end
%             end

        end
        
        function [pot , u , v, u_dx, u_dy, v_dx, v_dy] = induced_speed_on_many_points(VS, x, y)
            % Get Induced speed over multiple points
            
              pot  = zeros(size(x));
              u    = zeros(size(x));
              v    = zeros(size(x));
              u_dx = zeros(size(x));
              u_dy = zeros(size(x));
              v_dx = zeros(size(x));
              v_dy = zeros(size(x));

              for n_segment = 1:VS.n_segments
                  [pot_t , u_t , v_t, u_dx_t, u_dy_t, v_dx_t, v_dy_t] = ...
                    constant_strenght_vortex_segment_element(x , y, ...
                                                         VS.x_start(n_segment), VS.y_start(n_segment), ...
                                                         VS.x_end(n_segment)  , VS.y_end(n_segment), ...
                                                         VS.gamma(n_segment)); 
                   % Sum contribution of each segment on all points (we are
                   % summing arrays here, semi-efficient approach (we
                   % vectorize only one of two loops!)
                   pot  = pot  + pot_t;
                   u    = u    + u_t;
                   v    = v    + v_t;
                   u_dx = u_dx + u_dx_t;
                   u_dy = u_dy + u_dy_t;
                   v_dx = v_dx + v_dx_t;
                   v_dy = v_dy + v_dy_t;
              end
            
        end
    end
    
end

