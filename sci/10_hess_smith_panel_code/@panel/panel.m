classdef panel < handle
    %PANEL Class for storing and computing local information on a panel
    %   Code is and comented and self explanatory!
    %   Gael de Oliveira
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
        % Public properties
        x1      % First Point  x coordinate
        y1      % First Point  y coordinate
        
        x2      % Second Point x coordinate
        y2      % Second Point y coordinate
        
        u_norm_BC = 0 % Prescribed Normal Velocity
        
        gamma   % vortex strenght
        q       % source Strenght
        
        % Dependent properties
        
        xc      % Center Point x coordinate
        yc      % Center Point y coordinate
        
        l       % Panel Length
        
        sin_theta % Panel Inclination Angle sine
        cos_theta % Panel Inclination Angle cosine
        
        % Panel Base written in Original Base
        n_versor  % Normal base vector
        t_versor  % Tangential base vector               
        
        % Original Base written in Panel Base
        x_versor_local_base % x base versor
        y_versor_local_base % y base versor
        
        u_n_BC  % Prescribed Normal Velocity
        
        rel_exteriority = 1e-2 % Offset outside of panel for beta to be right
    end
    
    methods
        function p = panel(x1, y1 , x2, y2 , u_n_BC)
            % Class creator
            % Assign Prescribed fields
            p.x1 = x1;
            p.y1 = y1;
            
            p.x2 = x2;
            p.y2 = y2;
            
            % Set normal velocity boundary condition to specification or
            % zero if nothing specified
            if isempty(u_n_BC)
                p.u_n_BC = 0; 
            else
                p.u_n_BC = p.u_norm_BC;
            end
            
            % Initialize Strength of singularities
            p.gamma = 0;
            p.q = 0;        
            
            % Panel
            
            % Now compute everything needed
            update_panel(p)
        end
        
        function update_panel(p)
            % (re-)do Center Point coordinates (eg. if x1/2,y1/2 changed)
            make_geometry(p);
            
            % (re-)do                                             
        end
        
        function p = make_geometry(p)
        
        % Make Panel Length
        p.l  = sqrt((p.x1 - p.x2).^2 + (p.y1 - p.y2).^2);
                
        % Compute Panel Inclination Angles ( and store explicitly to avoid
        % singularities of atan)        
        p.sin_theta = (p.y2 - p.y1) / p.l;
        p.cos_theta = (p.x2 - p.x1) / p.l;
        
        % Define Original Base, before making Panel Base
        x_versor = [1 0];
        y_versor = [0 1];        
        % It is important to keep a right handed system (an error in the
        % tangential velocities seems to occur otherwise)
        p.n_versor =    p.sin_theta * x_versor - p.cos_theta * y_versor;
        p.t_versor =    p.cos_theta * x_versor + p.sin_theta * y_versor;
        
        
        % Now write original base in panel reference frame
        imat = inv([p.t_versor ; p.n_versor]);
        p.x_versor_local_base = imat(1,:);
        p.y_versor_local_base = imat(2,:);
        
        
        % Make Panel Midpoint and offset slightly outside for beta_ii to
        % ref
%         p.xc = (p.x1 + p.x2) / 2;
%         p.yc = (p.y1 + p.y2) / 2;
        
        offset = p.rel_exteriority * p.l;
        p.xc = (p.x1 + p.x2) / 2 + offset * p.n_versor(1);
        p.yc = (p.y1 + p.y2) / 2 + offset * p.n_versor(2);
                    
        end
        
        function r1 = r1_function(p , x , y)
            % Function returning the distance of point (x,y) in general
            % coordinates to point (x1, y1) of panel
            r1 = sqrt((x - p.x1).^2 + (y - p.y1).^2);                        
        end
        
        function r2 = r2_function(p , x , y)
            % Function returning the distance of point (x,y) in general
            % coordinates to point (x2, y2) of panel
            r2 = sqrt((x - p.x2).^2 + (y - p.y2).^2);                        
        end
        
        function beta = beta_function(p, x , y)
            % Function returns the angle beta for an arbitrary point (x,y)            
            % Using ATAN2(Y,X) to compute four quadrant arctangent
            % IMPLEMENT : When panel control point is requested supply Bij = pi as we
            % are interested in exterior problem
            
            % Compute distance vector to first node in global ref frame
            r1x = x - p.x1;
            r1y = y - p.y1;
            % Now move to local frame of reference
            [r1t r1n] = rewrite_global_vector_in_local_frame(p, r1x, r1y);
            % And compute auxilliary angle!
            v0 = atan2(r1n , r1t);
            
            r2x = x - p.x2;
            r2y = y - p.y2;
            % Now move to local frame of reference
            [r2t r2n] = rewrite_global_vector_in_local_frame(p, r2x, r2y);
            % And compute auxilliary angle!
            vl = atan2(r2n , r2t);
            
            beta = vl - v0;

        end
                
        function [u_x u_y] = rewrite_local_vector_in_global_frame(p, u_t, u_n)
            % Create global vector by linear combination of local base
            % elements
            u_global = u_t * p.t_versor + u_n * p.n_versor;
            
            % Extract elements and return!
            u_x = u_global(1);
            u_y = u_global(2);
        end
        
        function [u_t u_n] = rewrite_global_vector_in_local_frame(p, u_x, u_y)
            % Create global vector by linear combination of local base
            % elements            
            u_local  = u_x * p.x_versor_local_base  + u_y * p.y_versor_local_base;
            
            % Extract elements and return!
            u_t = u_local(1);
            u_n = u_local(2);
        end
        
        function [u_t u_n] = unit_source_induced_speeds_local_frame(p, x ,y)
            % Compute Distance to each panel extremity
            r1 = r1_function(p , x , y);
            r2 = r2_function(p , x , y);
            
            % Compute Beta Angle
            beta12 = beta_function(p, x , y);
            
            % Local speed components
            u_t = - 1/ (2*pi) * log(r2 ./ r1);
            u_n = beta12 / (2*pi());
        end
        
        function [u_x u_y] = unit_source_induced_speeds_global_frame(p, x ,y)
            % Get Induced Speed Vector in Panel frame of reference
            [u_t u_n] = p.unit_source_induced_speeds_local_frame(x ,y);
            
            % Convert to Global frame of reference
            [u_x u_y] = p.rewrite_local_vector_in_global_frame(u_t, u_n);
        end
                
        function [u_t u_n] = unit_vortex_induced_speeds_local_frame(p, x ,y)
            % Compute Distance to each panel extremity
            r1 = r1_function(p , x , y);
            r2 = r2_function(p , x , y);
            
            % Compute Beta Angle
            beta12 = beta_function(p, x , y);
            
            % Local speed components
            u_t = beta12 / (2*pi);
            u_n = 1/ (2*pi) * log(r2 ./ r1);
        end        
        
        function [u_x u_y] = unit_vortex_induced_speeds_global_frame(p, x ,y)
            % Get Induced Speed Vector in Panel frame of reference
            [u_t u_n] = unit_vortex_induced_speeds_local_frame(p, x ,y);
            
            % Convert to Global frame of reference
            [u_x u_y] = p.rewrite_local_vector_in_global_frame(u_t, u_n);                        
        end
        
    end    
end

