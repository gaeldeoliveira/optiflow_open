classdef constant_strenght_vortex_segment_2d_parallel < handle
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
        % Independent fields
        x_start             % Vortex segment start coordinates
        y_start             % 
        
        x_end               % Vortex segment end coodinates
        y_end               %
        
        gamma               % Vortex segment strenght
        
        % Independent fields (additional vortex pair)
        VP                  % Hook to vortex Pair Object
        
        % Dependent fields
        n_segments          % Number of segments
        
        x_center            % Vortex segment center coordinates
        y_center            %
        
        l                   % Vortex segment lenght
        
        x_l_unit_vector     % Directed Axis Unit vector components
        y_l_unit_vector     %
        
        x_n_unit_vector     % Normal Unit vector components
        y_n_unit_vector     % 
        
        N_threads = 4       % Number of threads to use when computing induced speed on many points
    end
    
    methods
        function VS = constant_strenght_vortex_segment_2d_parallel(x_start, y_start, x_end, y_end, gamma)
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
        end
        
        function [pot , u , v] = influence_on_single_point(VS, x, y)
            % Compute induced speed on single point
            [pot , u , v] = ...
                constant_strenght_vortex_segment_element(x , y, ...
                                                         VS.x_start, VS.y_start, ...
                                                         VS.x_end  , VS.y_end, ...
                                                         VS.gamma); 
        end
        
        function [pot , u , v] = induced_speed_on_single_point(VS, x, y)
            % Compute induced speed on single point
            
            % Get Single Point Influence
            [pot , u , v] = VS.influence_on_single_point(x, y);
            
            % Sum it up
            pot = sum(pot);
            u   = sum(u);
            v   = sum(v); 
            
            % Include contribution of lifting vortex pair, if applicable
            if not(isempty(VS.VP))
                % Compute contribution
                [pot_v, u_v , v_v] = VS.VP.induced_speed_on_many_points(x, y);
                % Add it
                pot = pot + pot_v;
                u   = u   + u_v;
                v   = v   + v_v;
            end
            
        end
        
        function [pot , u , v] = influence_on_many_points(VS, x, y) %#ok<STOUT,INUSD>
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
        
        function [pot , u , v] = induced_speed_on_many_points(VS, x, y)
            % Get Induced speed over multiple points
            
              % Allocate Return Arrays
              pot = zeros(size(x));
              u   = zeros(size(x));
              v   = zeros(size(x));
              
              % Allocate Temporary Thread Data Holders
              pot_cell = cell(1,VS.N_threads); pot_t_cell = cell(1,VS.N_threads);
              u_cell   = cell(1,VS.N_threads); u_t_cell   = cell(1,VS.N_threads);
              v_cell   = cell(1,VS.N_threads); v_t_cell   = cell(1,VS.N_threads);
              
              % Allocate thread data in data holder
              for n_thread = 1:VS.N_threads
                  pot_cell{n_thread} = pot;
                  u_cell{  n_thread} = u;
                  v_cell{  n_thread} = v;                  
              end
              
              % Spread Domain into threads
              thread_upper_boundaries = ceil((1:VS.N_threads) * VS.n_segments / VS.N_threads);
              thread_lower_boundaries = [1 (thread_upper_boundaries(1:(VS.N_threads-1))-1)];

              % Compute Contributions, which each thread handling a part of
              % the domain
              parfor n_thread = VS.N_threads
                  for n_segment = thread_lower_boundaries(n_thread):thread_upper_boundaries(n_thread)
                    [pot_t_cell{n_thread} , u_t_cell{n_thread} , v_t_cell{n_thread}] = ...
                                          constant_strenght_vortex_segment_element(x , y, ...
                                          VS.x_start(n_segment), VS.y_start(n_segment), ...
                                          VS.x_end(n_segment)  , VS.y_end(n_segment), ...
                                          VS.gamma(n_segment)); %#ok<PFOUS>
                    pot_cell{n_thread} = pot_cell{n_thread} + pot_t_cell{n_thread};
                    u_cell{  n_thread} = u_cell{  n_thread} + u_t_cell{  n_thread};
                    u_cell{  n_thread} = u_cell{  n_thread} + v_t_cell{  n_thread};
                  
                  end
              end
              
              % Now Sum up thread contributions
              for n_thread = VS.N_threads
                  pot = pot + pot_cell{n_thread};
                  u   = u   + u_cell{  n_thread};
                  v   = v   + v_cell{  n_thread};
              end
              
              % Include contribution of lifting vortex pair, if applicable
              if not(isempty(VS.VP))
                  % Compute contribution
                  [pot_v, u_v , v_v] = VS.VP.induced_speed_on_many_points(x, y);
                  % Add it 
                  pot = pot + pot_v;
                  u   = u   + u_v;
                  v   = v   + v_v;
              end
            
        end
    end
    
end

