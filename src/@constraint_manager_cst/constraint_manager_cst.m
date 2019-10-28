classdef constraint_manager_cst < handle
    % CONSTRAINT_MANAGER_CST class for handling constraints as accepted by
    % optimizers , ga_multiob, ga, and fmincon()
    %
    %   Detailed explanation goes here
    
    properties
        name
        id
        
        nvars
        
        lb0
        ub0
        
        lb_ext
        ub_ext
        
        A_ineq
        b_ineq
        
        A_eq
        b_eq
        
        extension_factor_ub_LE = 1;            % Greater value permits thicker airfoils (leading edge area)
        extension_factor_ub_M = 0.8;           % Greater value permits thicker airfoils (middle area)
        extension_factor_ub_TE = 1.5;          % Greater value permits thicker airfoils (trailing edge area)
        
        extension_factor_lb_LE = 1;            % Smaller value permits thinner airfoils
        extension_factor_lb_M = 1;             % Smaller value permits thinner airfoils
        extension_factor_lb_TE = 1;            % Smaller value permits thinner airfoils
        
        le_disc_constraint          = true;
        le_curv_max_diff            = 0.03;      % Maximum curvature discontinuity at LE between upper and lower side
        te_angle_constraint         = false;
        te_angle_max                = 8;         % Minimum trailing edge angle
        max_tc_constraint           = false;
        max_tc                      = 0.3;       % maximum thickness of airfoil     
        te_discontinuity_constraint = true;      % Control wether te_discontinuity_constraint is enabled (true means enabled)
        
        fixed_te_thickness_constraint = false    % fixed_te_thickness_constraint
        fixed_te_thickness = 0.003
        building_height_constraint = false       % Impose building_height equality constraint
        t_max_building_height = 0.3;             % chordwise position of maximum building_height
        building_height = 0.3;                   % specify building_height
        
        free_suction_area_constraints = false    % true means constraint for free suction area are set
        min_suction_length  = 0.1;               % Minimum length of suction region
        max_suction_length  = 0.4;               % Maximum length of suction region


        x2_min_suction       = 0.5;              % Minimum x/c for end of suction
        x2_max_suction       = 0.9;              % Maximum x/c for end of suction
        
        min_thickness_distribution
        tx_min_thickness_distribution
        
        ref_pop_array_bounds                     % Array containing reference population used for generation of upper and lower bounds
        names_fitted_bounds
        ref_pop_array_thickness                  % Array containing reference population used for generation of minimum thickness distribution
        names_fitted_thickness
        
        SC                                       % Handle to system context object
        SD                                       % Handle to shape definition object
        SF                                       % Handle to shape fit object
        
    end
    
    methods
        function CM = constraint_manager_cst(name, SC, SD, SF)
            CM.name = name;
            CM.SC = SC;
            CM.SD = SD;
            CM.SF = SF;
            CM.id   = now;
            
            % CM.suggest_constraints()
        end
        function determine_minimum_thickness_distribution(CM)
            
            % Now fit list of feasible airfoils
            CM.SF.obtain_airfoil_list_from_folder(CM.SC.feasible_airfoil_subdir);
            x_list_feasible = CM.SF.x_list_fitted;
            names_feasible  = CM.SF.names_fitted;
            
            CM.ref_pop_array_thickness  = x_list_feasible;     % Store loaded array for later reference
            CM.names_fitted_thickness   = names_feasible;      % Store loaded array names for later reference
            
            % Initialize Variables for minimum thickness distribution determination
            tx                  = zeros(201 , size(x_list_feasible,1));
            tz                  = zeros(201 , size(x_list_feasible,1));
            tx_thickness_camber = zeros(101 , size(x_list_feasible,1));
            thickness           = zeros(101 , size(x_list_feasible,1));
            camber              = zeros(101 , size(x_list_feasible,1));
            
            % Now analyse
            
            for n= 1:size(x_list_feasible,1)
                % Generate coordinates in classical and thickness camber form for each
                % airfoil in feasibility study folder
                [tx(:,n) tz(:,n)] = CM.SD.generate_coordinates(200 , x_list_feasible(n,:));
                [tx_thickness_camber(:,n) thickness(:,n) camber(:,n)] = CM.SD.generate_thickness_camber_coordinates(200 , x_list_feasible(n,:));
            end
                        
            % Now determine minimum thickness distribution
            CM.min_thickness_distribution    = min(thickness, [], 2);
            CM.tx_min_thickness_distribution = tx_thickness_camber(:,1);
        end
        function suggest_constraints(CM)
            % First determine minimum thickness distribution
            determine_minimum_thickness_distribution(CM)
            
            % Now determine bounds
            % Fit airfoils in bound database
            x_list = CM.SF.obtain_airfoil_list_from_folder(CM.SC.airfoil_subdir);
            % And store for latter reference
            CM.ref_pop_array_bounds = CM.SF.x_list_fitted;
            CM.names_fitted_bounds = CM.SF.names_fitted;                   
            % And really determine bounds from fitted airfoil list
            CM.make_bounds_from_airfoil_list(x_list);
            CM.make_linear_ineq_constraints();
            CM.make_linear_eq_constraints();
            
            l_order = CM.SD.parametrization_handles.lower.order;
            u_order = CM.SD.parametrization_handles.upper.order;
            N_dummy_arguments = CM.SD.N_dummy_parameters;
            
            CM.nvars = l_order + u_order + 1 + N_dummy_arguments;
        end
        function make_bounds_from_airfoil_list(CM, x_list)
            % This function suggests a set of bounds based on the list of
            % airofils in x_list and the extension factors defined in this
            % object's properties
            
            % Extract sizes of parametrizations!
            l_order = CM.SD.parametrization_handles.lower.order;
            u_order = CM.SD.parametrization_handles.upper.order;
            N_dummy_arguments = CM.SD.N_dummy_parameters;
            
            % This vector is used to limit expansion of bounds to cst shape
            % parameters, excluding trailind and dummy (non-shape)
            % parameters from this process
            
            %             % Original version allowing only single extension factor
            %             % per bound type
            %             extension_vector_base = [ones(1 , u_order + l_order) , ...
            %                 zeros(1,1) , zeros(1 , N_dummy_arguments)];
            %
            %             extension_factor_lb_vector = CM.extension_factor_lb * extension_vector_base;
            %             extension_factor_ub_vector = CM.extension_factor_ub * extension_vector_base;
            
            % Make variable expansion factor vector! The three user
            % specified values are fitted by a piecewise cubic polynomial
            x_ref = [0 0.5 1];
            y_ref_lb = [CM.extension_factor_lb_LE  CM.extension_factor_lb_M CM.extension_factor_lb_TE ];
            y_ref_ub = [CM.extension_factor_ub_LE  CM.extension_factor_ub_M CM.extension_factor_ub_TE ];
            
            x_interp_u_order = 0:1/(u_order-1):1;
            x_interp_l_order = 0:1/(l_order-1):1;
            x_interp = [x_interp_u_order , x_interp_l_order];
            
            extension_factor_lb_vector = [ interp1(x_ref , y_ref_lb , x_interp, 'pchip') , zeros(1,1) , zeros(1 , N_dummy_arguments)];
            extension_factor_ub_vector = [ interp1(x_ref , y_ref_ub , x_interp, 'pchip') , zeros(1,1) , zeros(1 , N_dummy_arguments)];
            
            % Turn list into array
            x_array = cell2mat(x_list);
            % Identify minimum values for airfoils in array
            CM.lb0    = min(x_array);
            % Identify maximum values for airfoils in array
            CM.ub0    = max(x_array);
            % Identify mean values for each coordinates over reference airfoil
            % population
            mean_x = mean(x_array);
            
            % Contruct extended airfoil bounds
            CM.lb_ext     =  mean_x - (mean_x - CM.lb0) .* extension_factor_lb_vector;
            CM.ub_ext     =  mean_x + (CM.ub0 - mean_x) .* extension_factor_ub_vector;
            
            % Add suction constraints if applicable
            if CM.free_suction_area_constraints == true
                % Suction coordinates defined as
                CM.lb_ext(end-1 : end)  = [CM.min_suction_length  , CM.x2_min_suction];
                CM.lb0(end-1 : end)     = [CM.min_suction_length  , CM.x2_min_suction];
                
                CM.ub_ext(end-1 : end)  = [CM.max_suction_length  , CM.x2_max_suction];
                CM.ub0(end-1 : end)     = [CM.max_suction_length  , CM.x2_max_suction];
            end
            
            % Always set lb of TE thickness to 0
            CM.lb_ext(l_order + u_order + 1) = 0;
            CM.lb0(l_order + u_order + 1)   = 0;
            
        end
        function make_linear_ineq_constraints(CM)
            %             le_curv_max_diff = CM.le_curv_max_diff;
            %             te_discontinuity_constraint = CM.te_discontinuity_constraint;
            %             max_suction_length = CM.max_suction_length;
            %             free_suction_area_constraints = CM.free_suction_area_constraints;
            % A version supporting multiple side sizes and constrained dummy arguments
            l_order = CM.SD.parametrization_handles.lower.order;
            u_order = CM.SD.parametrization_handles.upper.order;
            max_order = max(l_order , u_order);
            min_order = min(l_order , u_order);
            
            N_dummy_arguments = CM.SD.N_dummy_parameters;
            side_block_width = 1 + N_dummy_arguments;                   % 1 (te) + N_dummy_arguments
            
            % Block used to forbid crossing of airfoil sides (approximate linear constraint)      
            % Negative because constraint is of the type A*x < b and we want
            A_no_crossing_block = - [eye(max_order , u_order) , eye(max_order , l_order) , zeros(max_order, side_block_width)];
            b_no_crossing_block = zeros(max_order, 1);
                    
            % Block for controlling Leading Edge Discontinuity
            if CM.le_disc_constraint 
                A_le_disc_block = [[1 ,zeros(1 ,u_order-1) , -1 , zeros(1,l_order-1), zeros(1, side_block_width) ]; ...
                    [-1 ,zeros(1 ,u_order-1) , 1 , zeros(1,l_order-1), zeros(1, side_block_width)  ]] * CM.te_discontinuity_constraint;
                b_le_block      = CM.le_curv_max_diff * ones(2, 1);
            else
                A_le_disc_block = 1:0;
                b_le_block      = 1:0;
            end
            
%             % Block for controlling suction constraints (max length of suction region and min length)
%             %                     u_order + l_order + 1 (te)
%             suction_zeros_width = u_order + l_order + 1;
%             suction_order_line  = [zeros(1, suction_zeros_width) ,  1 , -1];
%             suction_length_line = [zeros(1, suction_zeros_width) , -1 ,  1];
%             if CM.free_suction_area_constraints == true
%                 A_suction_block       = [ suction_order_line ; suction_length_line ];
%                 b_suction_block       = [ 0 ; CM.max_suction_length ];
%             else
%                 A_suction_block       = 1:0;  % Make empty constraint block for suction as suction location constraints are disabled
%                 b_suction_block       = 1:0;
%             end
            
            % Block for minimum trailing edge angle constraint
            if CM.te_angle_constraint 
                A_te_angle_block    = [zeros(1 ,u_order-1),-1,zeros(1,l_order-1),-1,2,zeros(1, N_dummy_arguments)];
                b_te_angle_block    = -CM.te_angle_max*pi/180;
            else
                A_te_angle_block    = 1:0;
                b_te_angle_block    = 1:0;
            end
            
            if CM.max_tc_constraint 
                A_max_tc = [CM.make_building_height_line];
                b_max_tc = [CM.max_tc]; 
            else
                A_max_tc = 1:0;
                b_max_tc = 1:0;
            end
            
            % Build full A matrix and b vector
            CM.A_ineq = [ A_no_crossing_block ; A_le_disc_block ; A_te_angle_block ; A_max_tc];
            CM.b_ineq = [ b_no_crossing_block ; b_le_block ; b_te_angle_block ; b_max_tc];
            
        end
        function make_linear_eq_constraints(CM)
            % Function to make linear equality constraints. Make empty
            % constraints for now, but will do something more sophisticated
            % later, eventually allowing fixed point and value of maximum
            % thickness (check if doable in due time)
            
            % Initialize Constraints with no constraints, and add lines
            % according to user wishes!
            CM.A_eq = [];
            CM.b_eq = [];
            
            if CM.building_height_constraint == true
                CM.A_eq = [CM.A_eq ; CM.make_building_height_line() ; ...
                    CM.make_building_height_derivative_line()];
                CM.b_eq = [CM.b_eq ; CM.building_height ; 0];
                CM.fixed_te_thickness_constraint = true;
            end
            
            if CM.fixed_te_thickness_constraint == true
                N_dummy_parameters = CM.SD.N_dummy_parameters;
                Nu = CM.SD.parametrization_handles.upper.order;
                Nl = CM.SD.parametrization_handles.lower.order;
                
                fixed_te_thickness_line = ...
                    [ zeros( 1 ,Nu + Nl) , 1 , zeros( 1 ,N_dummy_parameters)];
                CM.A_eq = [CM.A_eq ; fixed_te_thickness_line];
                
                CM.b_eq = [CM.b_eq ; CM.fixed_te_thickness];
            end
            %            CM.A_eq = [CM.make_building_height_line()];
            %            CM.b_eq = [CM.building_height];
            
        end
        function [Aeq_line] = make_building_height_line(CM)
            
            % Extract necessary information
            Nl = CM.SD.parametrization_handles.lower.order;
            Nu = CM.SD.parametrization_handles.upper.order;
            
            M1u = CM.SD.parametrization_handles.upper.M1;
            M1l = CM.SD.parametrization_handles.lower.M1;
            
            N_dummy_parameters = CM.SD.N_dummy_parameters;
            
            t_p = CM.t_max_building_height;
            
            % Calculate class function
            C_t_p = sqrt(t_p) * (1 - t_p);
            
            % Build Canonic Polynomial Base Vector for upper part
            tp_u_vector = zeros(Nu, 1);
            for n_tp = 1:Nu
                tp_u_vector(n_tp) = t_p^(Nu-n_tp);              % Max degree = Nu -1 , and min degree equal 0!
            end
            
            % Get part of Aeq line corresponding to upper part
            Aeq_u_t = transpose(M1u * tp_u_vector);
            
            % Now for lower side
            % Build Canonic Polynomial Base Vector for lower part
            tp_l_vector = zeros(Nl, 1);
            for n_tp = 1:Nl
                tp_l_vector(n_tp) = t_p^(Nl-n_tp);              % Max degree = Nu -1 , and min degree equal 0!
            end
            
            % Get part of Aeq line corresponding to upper part
            Aeq_l_t = transpose(M1l * tp_l_vector);
            
            % Calculate trailing edge coefficient
            Aeq_te = t_p / C_t_p;
            
            % Now join the two parts together, and add slack for simulation parameters
            Aeq_line = C_t_p * [Aeq_u_t , Aeq_l_t, Aeq_te  , zeros(1, N_dummy_parameters)];
        end
        function [dAeq_line] = make_building_height_derivative_line(CM)
            
            % Extract necessary information
            Nu = CM.SD.parametrization_handles.upper.order;
            Nl = CM.SD.parametrization_handles.lower.order;
            
            M1u = CM.SD.parametrization_handles.upper.M1;
            M1l = CM.SD.parametrization_handles.lower.M1;
            
            Du = CM.SD.parametrization_handles.upper.D;
            Dl = CM.SD.parametrization_handles.lower.D;
            
            N_dummy_parameters = CM.SD.N_dummy_parameters;
            
            t_p = CM.t_max_building_height;
            
            % Calculate class function and its derivative
            C_t_p = sqrt(t_p) * (1 - t_p);
            dC_t_p_dt = 0.5 * (1 - t_p) / sqrt(t_p) - sqrt(t_p);
            
            % Make equation for upper side
            Dtu = transpose(Du);
            
            % Build Canonic Polynomial Base Vector for upper part
            tp_u_vector = zeros(Nu, 1);
            for n_tp = 1:Nu
                tp_u_vector(n_tp) = t_p^(Nu-n_tp);              % Max degree = Nu -1 , and min degree equal 0!
            end
            
            dSu_dt_line = transpose(M1u * Dtu * tp_u_vector);
            Su_t_line   = transpose(M1u * tp_u_vector);
            
            dSuC_dt_line = dC_t_p_dt * Su_t_line +  C_t_p * dSu_dt_line;
            
            % Now for lower Side
            Dtl = transpose(Dl);
            
            % Build Canonic Polynomial Base Vector for lower part
            tp_l_vector = zeros(Nl, 1);
            for n_tp = 1:Nl
                tp_l_vector(n_tp) = t_p^(Nl-n_tp);              % Max degree = Nl -1 , and min degree equal 0!
            end
            
            dSl_dt_line = transpose(M1l * Dtl * tp_l_vector);
            Sl_t_line   = transpose(M1l * tp_l_vector);
            
            dSlC_dt_line = dC_t_p_dt * Sl_t_line +  C_t_p * dSl_dt_line;
            
            % Now trailing edge coefficient
            dAeq_te = 1;
            
            % Now join everything together
            dAeq_line = [dSuC_dt_line , dSlC_dt_line, dAeq_te  , zeros(1, N_dummy_parameters)];
            
        end
    end
    
end


