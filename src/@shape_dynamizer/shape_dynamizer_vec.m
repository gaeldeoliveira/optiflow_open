classdef shape_dynamizer < handle
    %SHAPE_DYNAMIZER is a class used to support the simulation worker in
    % modifying the shape of the current sample. 
    % It can be used to represent:
    %       1. Flap Deflections with a smoothed 
    %       2. Apply a shape transformation to model VAWT rotational
    %       effects by considering the effect of unperturbed streamline
    %       curvature
    %       3. Refine the airfoil surface discretization in specific areas,
    %       for example to acccelerate convergence 
    %       4. Section Shape Deformation
    %       5. Any other thing you like!
    % The model may be refined further!
    
    properties
        current_trim_state      % Trim state vector        
        current_free_variables
        
        dynamic_variable_index  % Cell vector with list of field names that are variable
        fixed_variable_index    % Cell vector with list of field names that are fixed
        
        x_hinge = 0.8;          % Chordwise Stance Coordinate (x/c) of flap hinge
        z_hinge = 0;            % Normal Stance Coordinate (z/c) of flap hinge
        flap_knee  = 0.1;       % Smoothing Factor around Flap 
        
        flap_angle = 0;         % Flap angle in Degrees
        m_flap = 0;             % Flap hingeline slope to unflapped chordline
        
        SD                      % Shape Definition Object
        
    end
    
    methods
        function SDD = shape_dynamizer()
            
        end

        function generate_coordinates(SDD, n_points, parameters)
           % generate_coordinates method of shape_dynamizer, is an
           % overloaded version of the generate coordinates version of the
           % original shape_definition_cst class version
           
           
            
        end
        
        function h1 = h1_curvilinear_transformation(SDD, x, z)
           % Returns the h1 curvilinear transformation function      
           h1 = 1 - SDD.hingeline_curvature(x).*z;
        end
        
        function h2 = h2_curvilinear_transformation(SDD, x, z)
           % Returns the h2 curvilinear transformation function  
            h2 = 1;            
        end
                
        function s = T1_curvilinear_transformation(SDD, x, z)
           % Returns first component of T transformation
            s = x - SDD.hingeline_curvature_primitive(x).*z;            
        end
        
        function n = T2_curvilinear_transformation(SDD, x, z)
            % Returns second component of T transformation
            n = z + SDD.hingeline_deformation(x);
        end
        
        function [s, n] = T_curvilinear_transformation(SDD, x, z)
           % Applies to Transformation to return (s,n) when provided (x,z)
           s = T1_curvilinear_transformation(SDD, x, z - SDD.z_hinge);
           n = T2_curvilinear_transformation(SDD, x, z - SDD.z_hinge) + SDD.z_hinge;
        end
        
        function f = hingeline_deformation(SDD, x)
            % Returns the hingeline deformation at position x/c
            % Update Slope from Flap Angle
            m = SDD.flap_slope();
            delta = SDD.flap_knee;
            
            % Get Bounds of central curved region
            x_m_bound = SDD.x_hinge - SDD.flap_knee;
            x_p_bound = SDD.x_hinge + SDD.flap_knee;
           
            % Make conditions to process each streak separately
            c_m = (x<=x_m_bound);
            c_c = (and(x>=x_m_bound , x<x_p_bound));
            c_p = (x>=x_p_bound);            
            
            % Allocate space for return
            f = zeros(size(x));
            
            %%%% Take care of first streak (z_hinge = 0 for now! will be updated later)
            f(c_m) = x(c_m) * 0;            
            
            %%%% Take care of central curved streak
            % Get eta variable 
            eta = zeros(size(x));
            eta(c_c) = SDD.eta_function(x(c_c), delta);
            % Use p function (scaled g function) for curved region
            f(c_c) = SDD.p_function(eta(c_c));
            
            %%%% Take care of final straigh streak (curvature zero!)
            f(c_p) = (x(c_p) - SDD.x_hinge) * m; 
        end
        
        function k_f = hingeline_curvature(SDD, x)
            % Returns deformed hingeline curvature             
            
            % Update Slope from Flap Angle
            m = SDD.flap_slope();
            delta = SDD.flap_knee;
            
            % Get Bounds of central curved region
            x_m_bound = SDD.x_hinge - SDD.flap_knee;
            x_p_bound = SDD.x_hinge + SDD.flap_knee;
           
            x_m = x(x<=x_m_bound);
            x_c = x(and(x>=x_m_bound , x<x_p_bound));
            x_p = x(x>=x_p_bound);            
            
            %%%% Take care of first zero streak
            k_f_m = zeros(size(x_m));
            
            %%%% Take care of central curved streak
            % Get eta variable 
            eta = SDD.eta_function(x_c, delta);
                           
            % Compute hingeline_curvature_primitive for curved region
            k_f_c = (1./2*delta) .* (m*(- 6*eta.^2 + 6*eta)) ./ ...
                (m^2*(3*eta.^2 - 2*eta.^3).^2 + 1).^(3/2);
            
            
            %%%% Take care of final straigh streak (curvature zero!)
            k_f_p = zeros(size(x_p));
            
            
            %%%% Join all three streaks together
            k_f = [k_f_m , k_f_c , k_f_p];    
            
            
        end
        
        function k_f_prim = hingeline_curvature_primitive(SDD, x)
            % Returns deformed hingeline curvature primitive in x            
            
            % Update Slope from Flap Angle
            m = SDD.flap_slope();
            delta = SDD.flap_knee;
            
            % Get Bounds of central curved region
            x_m_bound = SDD.x_hinge - SDD.flap_knee;
            x_p_bound = SDD.x_hinge + SDD.flap_knee;
           
            x_m = x(x<=x_m_bound);
            x_c = x(and(x>=x_m_bound , x<x_p_bound));
            x_p = x(x>=x_p_bound);
            
            %%%% Take care of first zero streak
            k_f_prim_m = zeros(size(x_m));
            
            %%%% Take care of central curved streak
            % Get eta variable 
            eta = SDD.eta_function(x_c, delta);
                           
            % Compute hingeline_curvature_primitive for curved region
            k_f_prim_c = (m*(m^2*(- 2*eta.^3 + 3*eta.^2).^2 + 1).^(1/2) .* ... 
                    (- 2*eta.^3 + 3*eta.^2)) ./ ...
                        (4*eta.^6*m^2 - 12*eta.^5*m^2 + 9*eta.^4*m^2 + 1);
                    
                    
            %%%% Take care of final straigh streak                    
            eta = 1;
            k_f_prim_p_value = (m*(m^2*(- 2*eta.^3 + 3*eta.^2).^2 + 1).^(1/2)* ... 
                    (- 2*eta.^3 + 3*eta.^2)) ./ ...
                        (4*eta.^6*m^2 - 12*eta.^5*m^2 + 9*eta.^4*m^2 + 1);
                    
            k_f_prim_p = ones(size(x_p))*k_f_prim_p_value;
            
            
            %%%% Join all three streaks together
            k_f_prim = [k_f_prim_m , k_f_prim_c , k_f_prim_p];                    
        end
        
        function eta = eta_function(SDD, x, delta)
            % Returns the eta scaled coordinate from the x chordwise
            % coordinate, to generate curved region around flap hinge
            eta = (x - SDD.x_hinge) ./ (2 * delta) + 1/2;
        end
        
        function p = p_function(SDD, eta)
            % Returns the p function, unscaled version of the g function
            delta = SDD.flap_knee;
            m = SDD.flap_slope();
            p = 2 * delta * m * SDD.g_function(eta);
        
        end
        
        function g = g_function(SDD, eta)
            % Retrins the scaled g function determined from the BCs
            g = - 0.5 * eta.^4 + eta.^3;            
        end
        
        function m = flap_slope(SDD)
           m = tan(SDD.flap_angle * pi() / 180);
        end
    end
    
end

