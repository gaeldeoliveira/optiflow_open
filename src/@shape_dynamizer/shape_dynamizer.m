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
        dynamic_variable_list = {}              % Cell vector with list of field names that are variable
        dynamic_variable_index = []             % Array with index matching fieldnames with dummy_parameters elements (obtained from SD.breakdown_parameters(parameters))
        dynamic_variable_processors = {}        % Cell vector with processing functions setting how dummy_parameters are treated by each particular SDD object

        % % Mode Choice
        mode = 'flap';                          % If mode = 'flap', apply a flap transformation, if mode = 'thickness' , make constant thickness transformations
        
        % % Flap Parameters 
        x_hinge = 0.8;                          % Chordwise Stance Coordinate (x/c) of flap hinge
        z_hinge = 0;                            % Normal Stance Coordinate (z/c) of flap hinge
        flap_knee  = 0.1;                       % Smoothing Factor around Flap 
        
        flap_angle = 0;                         % Flap angle in Degrees
        m_flap = 0;                             % Flap hingeline slope to unflapped chordline
        
        % % Thickness Parameters 
        target_thickness = 0.25;                 % Target Thickness                
        
        
        SD                                      % Shape Definition Object
        relative_z_hinge = 0;                   % A scaled coordinate for z_hinge to lie inside airfoil in all cases
        flag_relative_z_hinge_dominance = true; % If true then z_hinge is determined from relative_z_hinge, otherwise the other way around        
    end
    
    methods
        function SDD = shape_dynamizer(varargin)
            % Constructs shape_dynamizer object and hooks it to SD object
            % if it supplied
            
            % Hook!
            if not(isempty(varargin))
                SDD.SD = varargin{1};
                SDD.SD.SDD = SDD;
            end
        end
        
        function [ts, tn] = dynamize_coordinates(SDD, parameters, tx, tz)
            % generate_coordinates method of shape_dynamizer, is an
            % overloaded version of the generate coordinates version of the
            % original shape_definition_cst class version
            
            % Start by parsing case information, inserting it into object
            % fields and making them consistent to define coordinate
            % transformation completely
            SDD.update_dynamic_fields(parameters);
            
            % Now make coordinate transformation
            if strcmp(SDD.mode, 'flap')
                [ts, tn] = SDD.T_curvilinear_transformation(tx, tz);
            end
            if strcmp(SDD.mode, 'thickness')
                [ts, tn] = SDD.thickness_projection_transformation(parameters, tx, tz);
            end
        end
        
        function update_dynamic_fields(SDD,parameters)
            % Parses airfoil/case definition parameters into SDD object
            % fields upon request during cost function evaluation
            
            % Check wether dynamic_variable_processors are
            % initialized/defined by user or external source
            if isempty(SDD.dynamic_variable_processors)
                % If not, initialize them
               SDD.initialize_dynamic_variable_processors();
            end
            
            % Extract dummy (non-CST-shape) parameters
            [parameters_upper , parameters_lower , z_te_upper , ...
                z_te_lower, dummy_parameters] = SDD.SD.breakdown_parameters(parameters);
            
            
            % Process/Parse dummy parameters into SDD fields
            for n = 1:length(SDD.dynamic_variable_list)
                SDD.(SDD.dynamic_variable_list{n}) = ...
                    SDD.dynamic_variable_processors{n}(...
                    dummy_parameters(SDD.dynamic_variable_index(n)));
            end
            
            % Update dependent fields to restore consistency now that input
            % fields are edited
            SDD.update_dependent_fields(parameters)
        end
        
        function update_dependent_fields(SDD, parameters)
            % Updates dependent fields to restore consistency after dynamic
            % fields are updated!
            
            % Get thickness and camber at x_hinge
            [thickness, camber] = SDD.SD.thickness_camber_at_tx(SDD.x_hinge, parameters);
            
            % Update one of the two dependent variables depending on flag_relative_z_hinge_dominance
            if SDD.flag_relative_z_hinge_dominance == true
                % If dominance flag on, then Update z_hinge from relative_z_hinge
                SDD.z_hinge = camber + 0.5 * thickness * SDD.relative_z_hinge;
            else
                % Otherwise Update relative_z_hinge from z_hinge
                SDD.relative_z_hinge = 2 * (SDD.z_hinge - camber) / thickness;
            end
            
        end
        
        function initialize_dynamic_variable_processors(SDD)
            % Initializes dynamic_variable_processors with unit functions
            % if undefined upon SDD.update_dynamic_fields() call
            
            % Allocate space
            SDD.dynamic_variable_processors = cell(size(SDD.dynamic_variable_list));
            
            % Not set each dynamic_variable_processor as a unit function!
            for n = 1:length(SDD.dynamic_variable_list)
                SDD.dynamic_variable_processors{n} = @(x) x;
            end
        end
        
        function [ts, tn] = thickness_projection_transformation(SDD, parameters, tx, tz)
            [tx_t, thickness, ~] = SDD.SD.generate_thickness_camber_coordinates(200, parameters);
            
            actual_thickness = max(thickness);
            
            scaling = SDD.target_thickness / actual_thickness;
            
            ts = tx;
            tn = tz * scaling;
            
        end
                
        function h1 = h1_curvilinear_transformation(SDD, x, z)
           % Returns the h1 curvilinear transformation function      
           h1 = 1 - SDD.hingeline_curvature(x).*z;
        end
        
        function h2 = h2_curvilinear_transformation(SDD, x, z) %#ok<INUSD>
           % Returns the h2 curvilinear transformation function  
            h2 = 1;            
        end
                
        function s = T1_curvilinear_transformation(SDD, x, z)
           % Returns first component of T transformation
           
           % Original Expression (consistent with derivation)
            % s = x - SDD.hingeline_curvature_primitive(x).*z;
            
            % New expression, I feel it is like this but somehow miss the
            % step where this starts making sense in the derivation. Maybe
            % the original transformation was conformal but not lenght
            % preserving, and here I need somehow both, but must choose to
            % a certain extent!
            s = x - SDD.hingeline_curvature_primitive(x).*z;
            % I understood it now, the normalization should be applied to
            % the hingeline deformation function itself, and the
            % derivatives and curvature that follow from there, which would
            % be completely length preserving. The above expression is an
            % approximation with vanishing error for small deflection
            % angles. It will stay like this for now, as it would give a
            % significant amount of work to rededuce the expressions and
            % adjust everything!
            
        end
        
        function n = T2_curvilinear_transformation(SDD, x, z)
            % Returns second component of T transformation
            n = z + SDD.hingeline_deformation(x);
        end
        
        function [s, n] = T_curvilinear_transformation(SDD, x, z)
           % Applies to Transformation to return (s,n) when provided (x,z)
           tx = x .* SDD.hingeline_normalization_element(x);
           s = T1_curvilinear_transformation(SDD, tx, z - SDD.z_hinge);
%            s = tx - SDD.hingeline_curvature_primitive(x).*z;          % %            Alternative solution that seems to make no difference
           s = tx - SDD.hingeline_curvature_primitive(tx).*z;
           n = T2_curvilinear_transformation(SDD, tx, z - SDD.z_hinge) + SDD.z_hinge;
           
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
        
        function norm1f = hingeline_normalization_element(SDD, x)
%             dnorm1f = 1 ./ sqrt(1 + (SDD.hingeline_deformation_x_derivative(x)).^2);
%              norm1f = zeros(size(x));
%              for n = 1:length(x)
%                  norm1f(n) = quadgk(@SDD.hingeline_normalization_element_x_derivative , 0, x(n)) ;                 
%              end
%              


             norm1f = ones(size(x));
             m = SDD.flap_slope();
             condition = (x-SDD.x_hinge)>0;
%              norm1f(not(condition)) = x(not(condition));
             norm1f(condition) = (SDD.x_hinge + (x(condition)-SDD.x_hinge) / sqrt(1+m^2)) ./ x(condition);
 %             % n1f = SDD.hingeline_deformation_x_derivative(x) ./ sqrt(1 + (SDD.hingeline_deformation_x_derivative(x)).^2);
%             norm1f = 0;
        end
        
        function dnorm1f = hingeline_normalization_element_x_derivative(SDD, x)
            dnorm1f = 1 ./ sqrt(1 + (SDD.hingeline_deformation_x_derivative(x)).^2);
        end
                
        function d1f = hingeline_deformation_x_derivative(SDD, x)
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
            d1f = zeros(size(x));
            
            %%%% Take care of first streak (d1f = 0)
            d1f(c_m) = x(c_m) * 0;            
            
            %%%% Take care of central curved streak
            % Get eta variable 
            eta = zeros(size(x));
            eta(c_c) = SDD.eta_function(x(c_c), delta);
            % Use p function (scaled g function) for curved region
            d1f(c_c) = SDD.d1p_function(eta(c_c));
            
            %%%% Take care of final straight streak (d1f = m)
            d1f(c_p) = c_p(c_p) * m; 
        end
        
        function k_f = hingeline_curvature(SDD, x)
            % Returns deformed hingeline curvature             
            
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
            k_f = zeros(size(x));
            
            %%%% Take care of first zero streak
            k_f(c_m) = x(c_m) * 0;
            
            %%%% Take care of central curved streak
            % Get eta variable 
            eta = zeros(size(X));
            eta(c_c) = SDD.eta_function(x(c_c), delta);
                           
            % Compute hingeline_curvature_primitive for curved region
            k_f(c_c) = (1./2*delta) .* (m*(- 6*eta(c_c).^2 + 6*eta(c_c))) ./ ...
                (m^2*(3*eta(c_c).^2 - 2*eta(c_c).^3).^2 + 1).^(3/2);
            
            %%%% Take care of final straigh streak (curvature zero!)
            k_f(c_p) = x(c_p) * 0;                        
        end
        
        function k_f_prim = hingeline_curvature_primitive(SDD, x)
            % Returns deformed hingeline curvature primitive in x            
            
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
            k_f_prim = zeros(size(x));
            
            %%%% Take care of first zero streak
            k_f_prim(c_m) = x(c_m)*0;
            
            %%%% Take care of central curved streak
            % Get eta variable 
            eta = zeros(size(x));
            eta(c_c) = SDD.eta_function(x(c_c), delta);
                           
            % Compute hingeline_curvature_primitive for curved region
            k_f_prim(c_c) = (m*(m^2*(- 2*eta(c_c).^3 + 3*eta(c_c).^2).^2 + 1).^(1/2) .* ... 
                    (- 2*eta(c_c).^3 + 3*eta(c_c).^2)) ./ ...
                        (4*eta(c_c).^6*m^2 - 12*eta(c_c).^5*m^2 + 9*eta(c_c).^4*m^2 + 1);
                    
                    
            %%%% Take care of final straigh streak                    
            eta = 1;
            k_f_prim_p_value = (m*(m^2*(- 2*eta.^3 + 3*eta.^2).^2 + 1).^(1/2)* ... 
                    (- 2*eta.^3 + 3*eta.^2)) ./ ...
                        (4*eta.^6*m^2 - 12*eta.^5*m^2 + 9*eta.^4*m^2 + 1);

            k_f_prim(c_p) = c_p(c_p) *k_f_prim_p_value;
            
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
        
        function d1p = d1p_function(SDD, eta)
            % Returns the x derivative of the p function at x = x(eta), unscaled version of the g function            
            m = SDD.flap_slope();
            d1p = m * SDD.d1g_function(eta);        
        end
        
        function g = g_function(SDD, eta)
            % Returns the scaled g function determined from the BCs
            g = - 0.5 * eta.^4 + eta.^3;            
        end
        
        function d1g = d1g_function(SDD, eta)
            % Returns the eta derivative of the scaled g function determined from the BCs
            d1g = - 2 * eta.^3 + 3 * eta.^2;
        end
        
        function m = flap_slope(SDD)
           m = tan( - SDD.flap_angle * pi() / 180);     % Slope is with negative flap angle to follow usual convention in aerodynamics (lift goes up when flap angle increases!)! 
        end
    end
    
end

