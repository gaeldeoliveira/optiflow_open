classdef airfoil_panel_set
    %AIRFOIL_PANEL_SET
    % Manager class (non-handle) for set of airfoil panels
    % Gael de Oliveira 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                     %
    %   ActiHS  :   A simple modular Panel Code based on a free           %
    %               interpretation of the Hess-Smith method in velocity   %
    %               components                                            %
    %                                                                     %
    %   Usage   :   Standalone with script, for sail optimization, within %
    %               the kirikou-dogoro actuator codes or other codes and  %
    %               derivatives from the author                           %
    %                                                                     %
    %   Date    :   April 2011 to March 2017                              %
    %   Author  :   Gael de Oliveira                                      %
    %                                                                     %
    %   License :   MIT, as the rest of this repository 		    %
    %                                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        panel_list
        
        N_panel
        N_nodes 
        
        q_vector
        Gamma_vector        
    end
    
    methods
        function aps = airfoil_panel_set(panel_list_input)
            aps.panel_list = panel_list_input;
            aps.N_panel = length(panel_list_input);
            aps.N_nodes = aps.N_panel + 1;
        end
        
        function [u_x_line, u_y_line] = source_induced_speeds_line(aps , x , y)
            u_x_line = zeros(1,aps.N_panel);
            u_y_line = zeros(1,aps.N_panel);
            
            for n_panel = 1:aps.N_panel
                [u_x, u_y] = aps.panel_list{n_panel}.unit_source_induced_speeds_global_frame(x,y);
                u_x_line(n_panel) = u_x;
                u_y_line(n_panel) = u_y;
            end
            
        end
        
        function [u_x_line, u_y_line] = vortex_induced_speeds_line(aps , x , y)
            u_x_line = zeros(1,aps.N_panel);
            u_y_line = zeros(1,aps.N_panel);
            
            for n_panel = 1:aps.N_panel
                [u_x, u_y] = aps.panel_list{n_panel}.unit_vortex_induced_speeds_global_frame(x,y);
                u_x_line(n_panel) = u_x;
                u_y_line(n_panel) = u_y;
            end            
        end
        
        function [u_x_M, u_y_M] = source_induced_speeds_matrices(aps)
            % Influence Matrix of sources on 
            u_x_M = zeros(aps.N_panel , aps.N_panel);
            u_y_M = zeros(aps.N_panel , aps.N_panel);
            % Each line corresponds to a different control point
            for n_c = 1:aps.N_panel
                % Identify control point of current panel
                xc = aps.panel_list{n_c}.xc;
                yc = aps.panel_list{n_c}.yc;
                
                % And generate source speed contributions at that point.
                [u_x_line, u_y_line] = source_induced_speeds_line(aps , xc , yc);
                
                u_x_M(n_c , :) = u_x_line;
                u_y_M(n_c , :) = u_y_line;
            end
            
        end
        
        function [u_x_M, u_y_M] = vortex_induced_speeds_matrices(aps)
            % Influence Matrix of sources on 
            u_x_M = zeros(aps.N_panel , aps.N_panel);
            u_y_M = zeros(aps.N_panel , aps.N_panel);
            % Each line corresponds to a different control point
            for n_c = 1:aps.N_panel
                % Identify control point of current panel
                xc = aps.panel_list{n_c}.xc;
                yc = aps.panel_list{n_c}.yc;
                
                % And generate source speed contributions at that point.
                [u_x_line, u_y_line] = vortex_induced_speeds_line(aps , xc , yc);
                
                u_x_M(n_c , :) = u_x_line;
                u_y_M(n_c , :) = u_y_line;
            end
        end
        
        function [sin_theta_M, cos_theta_M] = trigonometic_matrices(aps)
            % Make matrices for transformation of u_x u_y speeds into local
            % frame of reference
            % Diagonal matrices with sin_theta and cos_theta for each panel
            
            % Prealocate diagonal entry vectors
            sin_theta_vec = zeros(size(aps.panel_list));
            cos_theta_vec = zeros(size(aps.panel_list));
            
            % Now fill them
            for n_c = 1:aps.N_panel
                sin_theta_vec(n_c) = aps.panel_list{n_c}.sin_theta;
                cos_theta_vec(n_c) = aps.panel_list{n_c}.cos_theta;
            end
            
            % And generated diagonal matrices we were aiming for
            sin_theta_M = diag(sin_theta_vec);
            cos_theta_M = diag(cos_theta_vec);            
        end        
        
        function [u_n_M] = normal_speeds_matrix(aps , u_x_M , u_y_M)
            
            [sin_theta_M, cos_theta_M] = trigonometic_matrices(aps);            
            
            u_n_M = sin_theta_M * u_x_M - cos_theta_M * u_y_M;
        end
        
        function [u_t_M] = tangential_speeds_matrix(aps , u_x_M , u_y_M)
            
            [sin_theta_M, cos_theta_M] = trigonometic_matrices(aps);
            
            u_t_M = cos_theta_M * u_x_M + sin_theta_M * u_y_M;
        end
        
        function [u_n_BC_vector] = normal_BC_condition(aps)
            % Makes vertical vector with list of normal BC velocities
            
            % Allocate memory
            u_n_BC_vector = zeros(aps.N_panel,1);
            
            % Fill vector in
            for n_panel = 1:aps.N_panel
                u_n_BC = aps.panel_list{n_panel}.u_n_BC;
                u_n_BC_vector(n_panel) = u_n_BC;
            end
            % Done!
        end
        
        function [x_c_column, y_c_column] = centerpoint_locations_columns(aps)
            % Returns locations of centerpoints in two column vectors
            
            % Allocate memory
            x_c_column = zeros(aps.N_panel,1);
            y_c_column = zeros(aps.N_panel,1);
            
            % Now fill vectors with relevant data
            for n_panel = 1:aps.N_panel
                x_c_column(n_panel) = aps.panel_list{n_panel}.xc;
                y_c_column(n_panel) = aps.panel_list{n_panel}.yc;                
            end            
            % Done, return!
        end
        
        function [n_versor_x, n_versor_y] = normal_versors_column(aps)
           % Returns ordered list of panel normal versors 
           
           % Allocate memory
           n_versor_x = zeros(aps.N_panel,1);
           n_versor_y = zeros(aps.N_panel,1);
           
           % Move on to fill in versors
           for n_panel = 1:aps.N_panel
               n_versor = aps.panel_list{n_panel}.n_versor;
               n_versor_x(n_panel) = n_versor(1);
               n_versor_y(n_panel) = n_versor(2);
           end            
           % Done, return!
        end
        
        function [l_column] = panel_lengths_column(aps)
           % Returns ordered list of panel lengths
           
           % Allocate memory
           l_column = zeros(aps.N_panel,1);           
           
           % Move on to fill in lengths
           for n_panel = 1:aps.N_panel               
               l_column(n_panel) = aps.panel_list{n_panel}.l;
           end            
           % Done, return!
        end
        
    end
    
end


