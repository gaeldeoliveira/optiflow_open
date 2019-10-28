classdef shape_definition_constant_thickness < shape_definition_cst
    % SHAPE_DEFINITION This class stores, aquires and converts airfoil shape
    % definitions
     
    properties
        N_dummy_parameters  = 0;    % Number of parameters left at the end of the vector that are not related to the shape definition % a(end+1:end+N_dummy_parameters) = 0
        SDD = [];                   % Shape Dynamizer Object, Leave empty if unnecessary
    end
    
    properties(SetAccess = private) 
        
    end
    
    methods
        function SD = shape_definition_cst(parametrization_handle_upper,  parametrization_handle_lower, N_dummy_parameters, name)
            SD.parametrization_handles.upper = parametrization_handle_upper;
            SD.parametrization_handles.lower = parametrization_handle_lower;
            SD.N_dummy_parameters = N_dummy_parameters;
            SD.name = name;
            SD.id = now();    
        end
        
        function tz = tz_upper(SD, tx, parameters_upper , z_te_upper)
            C_t  = SD.parametrization_handles.upper.C_t;
            tz =   z_side(tx, parameters_upper, z_te_upper , C_t);
            % Mimi L'Ennui (Bobino 80) - Renaud
        end
        
        function tz = tz_lower(SD, tx, parameters_lower , z_te_lower)
            C_t  = SD.parametrization_handles.upper.C_t;            
            tz = - z_side(tx, parameters_lower, z_te_lower , C_t);
        end
        
        function [parameters_upper , parameters_lower , z_te_upper , z_te_lower , dummy_parameters] = breakdown_parameters(SD , parameters)
%             % Code for same order on upper and lower side
%             N = length(parameters) - SD.N_dummy_parameters;                       
%             parameters_upper = parameters(1:((N-1)/2));
%             parameters_lower = parameters(((N-1)/2+1):(N-1));
            
            % Code compatible with different order on upper and lower side
            u_order = SD.parametrization_handles.upper.order;
            l_order = SD.parametrization_handles.lower.order;
            parameters_upper = parameters(1:u_order);
            parameters_lower = parameters((u_order+1):(u_order+l_order));
                                    
            z_te = parameters(u_order+l_order+1);
            z_te_upper = z_te;
            z_te_lower = z_te;
            % Chanson pour Pierrot - Renaud
            dummy_parameters = parameters((u_order+l_order+2):end);
        end
        
        function [tx, tz] = generate_coordinates(SD, n_points, parameters)
            % This function should be greatly improved before it is
            % documented in its final version. 
            % Currently refinement only takes place at leading edge. It
            % would be however nice to have adaptive refinement where
            % separation occurs:
            %       TE on suction side
            %       Near laminar sep. bubble on pressure side
                                    
            % If number of points is not odd, add one to ensure full domain
            % coverage
            if floor(n_points/2) == n_points/2
                n_points = n_points + 1;
            end
            
            % Generate abcissa vector (number of points must be pair)
            tx_lower = 1-cos((0:1/((n_points-1)/2):1) * pi/2);
            tx_upper = fliplr(tx_lower(2:end));
            
            % Breakdown parameter vector
            [parameters_upper , parameters_lower , z_te_upper , z_te_lower] = SD.breakdown_parameters(parameters);
            
            % Generate Ordinate vectors
            tz_upper = SD.tz_upper(tx_upper , parameters_upper , z_te_upper);
            tz_lower = SD.tz_lower(tx_lower , parameters_lower , z_te_lower);
            
            % Concatenate Vectors to obtain airfoil coordinates!
            tx = transpose([tx_upper , tx_lower]);
            tz = transpose([tz_upper , tz_lower]);
            
            % Dynamize Coordinates if SDD object is present
            if not(isempty(SD.SDD))
                [tx, tz] = SD.SDD.dynamize_coordinates(parameters, tx, tz);
            end
        end
        
        function [tx, thickness, camber] = generate_thickness_camber_coordinates(SD, n_points, parameters)
            % This function should be greatly improved before it is
            % documented in its final version. 
            % Currently refinement only takes place at leading edge. It
            % would be however nice to have adaptive refinement where
            % separation occurs:
            %       TE on suction side
            %       Near laminar sep. bubble on pressure side
                                    
            % If number of points is not odd, add one to ensure full domain
            % coverage
            if floor(n_points/2) == n_points/2
                n_points = n_points + 1;
            end
            
            % Generate abcissa vector (number of points must be pair)
            tx = 1-cos((0:1/((n_points-1)/2):1) * pi/2);            
            
            % Breakdown parameter vector
            [parameters_upper , parameters_lower , z_te_upper , z_te_lower] = SD.breakdown_parameters(parameters);
            
            % Generate Ordinate vectors
            tz_upper = SD.tz_upper(tx, parameters_upper , z_te_upper);
            tz_lower = SD.tz_lower(tx, parameters_lower , z_te_lower);
            
            % Concatenate Vectors to obtain airfoil coordinates!
            thickness = tz_upper - tz_lower;
            camber    = (tz_upper + tz_lower) / 2;
        end
        
        function [thickness, camber] = thickness_camber_at_tx(SD, tx, parameters)
            % Returns thickness and camber at stance tx for airfoil
            % defined by parameters
            
            % Breakdown parameter vector
            [parameters_upper , parameters_lower , z_te_upper , z_te_lower] = SD.breakdown_parameters(parameters);
            
            % Generate Ordinate vectors
            tz_upper = SD.tz_upper(tx, parameters_upper , z_te_upper);
            tz_lower = SD.tz_lower(tx, parameters_lower , z_te_lower);
            
            % Concatenate Vectors to obtain airfoil coordinates!
            thickness = tz_upper - tz_lower;
            camber    = (tz_upper + tz_lower) / 2;
            
        end
        
        function airfoil_figure = plot_airfoil(SD, parameters, varargin)
            [tx, tz] = SD.generate_coordinates(200, parameters);
            
            if isempty(varargin)
                figure()
                airfoil_figure = figure();
            else
                airfoil_figure = varargin{1};
            end
            
            % Check if a color was specified
            if length(varargin) >= 2                
                plot_color = varargin{2};
            else
                % If none, set as black
                plot_color = [0 0 0]; 
            end
            
            % Now, finally: Plot!!!
            figure(airfoil_figure)
            hold on
            plot(tx, tz, 'Color', plot_color);
            hold off
            
            % Set equal scale axes and grid
            axis equal
            grid on
            
            % Set labels
            xlabel('x/c')
            ylabel('z/c')
            
            % Restons Amants - Maxime le Forestier
        end
    end
    
end

% Soan - Emily