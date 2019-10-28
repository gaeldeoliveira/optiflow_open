classdef inviscid_panel_case_multielement < handle
    %INVISCID_PANEL_CASE_MULTIELEMENT manager class for inviscid flow case
    %   with multiple lifting elements interacting with each other
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
        
        % Default Parameters
        u_inf = 1;              % Specify Flow base speed
        alpha = 0 * pi / 180;   % Angle of attack        
        c_over_R = 0.1;         % Now specify rotational effects
        x_cp = 1/4;             % Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
        
        rotation = false;       % Is effect of rotation accounted for ? true = yes
        
        % Airfoil Coordinates 
        px                      % x-coordinates, usually in x/c units, with R/Xfoil point ordering (TEup->LE->TEdown)
        py                      % y-coordinates, usually in y/c units, with R/Xfoil point ordering (TEup->LE->TEdown)
        
        te_first_indices        % Vector with list of first points on airfoils (coinciding with TE) / index given in matrices reference
        te_last_indices         % Vector with list of lasst points on airfoils (coinciding with TE) / index given in matrices reference
        
        te_first_indices_input  % Vector with list of first points on airfoils (coinciding with TE) / index given in inputs reference
        te_last_indices_input   % Vector with list of last points on airfoils (coinciding with TE)  / index given in inputs reference                
        
        N_elements              % Indicates number of elements on multielement airfoil
                
        % Dependent parameters
        R                       % (assume chord is equal to 1, LE at origin and TE at (1,0))        
        x_center                % Center of Rotation Location (x coordinate)
        y_center                % Center of Rotation Location (y coordinate)
        p_center                % In vector form
        
        L_versor
        
        % Handle Functions
        u_x_theta
        u_y_theta
        
        % Geometry Holders
        aps
        
        % Data exchange elements
        u_n_BC_vector
        u_v_t_M
        u_s_t_M 
        u_t_M
        u_n_M
        u_inf_x_vec             % x-component          of forcing (free-stream) speed at panel reference points
        u_inf_y_vec             % y-component          of forcing (free-stream) speed at panel reference points
        u_inf_n_vec             % normal-component     of forcing (free-stream) speed at panel reference points
        u_inf_t_vec             % tangential-component of forcing (free-stream) speed at panel reference points
        
        % System Matrices
        A                       % The big matrix! for a system in Ax = b form
        b                       % The loading column
        q_vector                % The solution vector
        
        % Answers (Cp, Lift and such things)
        L_Cp                    % Lift Computed from Pressure Coefficients
        L_gamma                 % Lift Computed from Circulation
        D_Cp                    % Drag Computed from Pressure Coefficients
        
        res                     % Primary Residual   : RMS deviation to impermeability condition
        res2                    % Secondary Residual : 
        
        u_t_vec                 % Tangential Velocity Vector
        cp                      % Vector of pressure coefficients at px_middle stances
        cp_plot                 % Vector of pressure coefficients shifted to 0 at free stream
        px_middle               % x stances of mid panel points 
        py_middle               % y stances of mid panel points 
        
        F_x                     % Force in X direction
        F_y                     % Force in Y direction
    end
    
    methods
        function ipc = inviscid_panel_case_multielement (px , py , te_first_indices_input , te_last_indices_input)
            ipc.px = px;
            ipc.py = py;
            ipc.te_first_indices_input = te_first_indices_input;
            ipc.te_last_indices_input  = te_last_indices_input;                        
            
            ipc.N_elements = length(te_first_indices_input);
            
            generate_A_matrix(ipc);
            
            [ipc.px_middle, ipc.py_middle] = ipc.aps.centerpoint_locations_columns();
        end
        
        function generate_solution(ipc)
            update_dependent_inputs(ipc)
            make_rotation_freestream_function_handles(ipc);           
            generate_b_column(ipc);
            solve(ipc);
            post_process(ipc);            
        end
        
        function update_dependent_inputs(ipc)
            ipc.R = 1 / ipc.c_over_R;

            % Generate Center of Rotation Location
            ipc.x_center = ipc.x_cp;
            ipc.y_center = -ipc.R;
            % And put into vector
            ipc.p_center = [ipc.x_center , ipc.y_center];
        end
        
        function make_rotation_freestream_function_handles(ipc)
            % Move on to define rotation
            % Make function for r vector from point (x,y) to center
            r_to_center = @(x,y) [x,y] - ipc.p_center;
            % Norm of r vector, distance between point (x,y) and center of rotation
            r_norm = @(x,y) norm(r_to_center(x,y));
            % Make azimuthal direction versor using a matrix with unit determinant (is
            % direction right?)
            theta_versor = @(x,y) r_to_center(x, y) * [0 -1 ; 1 0] / r_norm(x, y);
            % And go on to define speed vector;
            u_vec_theta = @(x,y) theta_versor(x,y) * ipc.u_inf * r_norm(x,y) / ipc.R;
            % Now provide for single coordinate extraction so that vectorizatin with
            % arrayfun is easy
            ipc.u_x_theta = @(x,y) u_vec_theta(x,y) * [1 ;0];
            ipc.u_y_theta = @(x,y) u_vec_theta(x,y) * [0 ;1];            
        end
        
        function generate_A_matrix(ipc)
            %% Now make panels
            % Original Version            
%             % Allocate memory
%             ipc.panel_list = cell(length(ipc.px)-1, 1);
%             % Now run loop to create every panel
%             for n_panel = 1:length(ipc.panel_list)
%                 ipc.panel_list{n_panel} = panel(ipc.px(n_panel) , ipc.py(n_panel) , ...
%                     ipc.px(n_panel+1) , ipc.py(n_panel+1) , []);
%             end
            % Multielement Version
            % Allocate memory
            panel_list = cell(length(ipc.px)-ipc.N_elements, 1);
            % Run loop to create every panel            
            n_panel_abs = 1;
            % For each element
            for n_element = 1:ipc.N_elements
                % Run through each panel, with local index, n_panel, and
                % stuff things into panel list with global index n_panel_abs
                ipc.te_first_indices(n_element) = n_panel_abs;
                for n_panel = ipc.te_first_indices_input(n_element):(ipc.te_last_indices_input(n_element)-1)
                    panel_list{n_panel_abs} = panel(ipc.px(n_panel) , ipc.py(n_panel) , ...
                        ipc.px(n_panel+1) , ipc.py(n_panel+1) , []);
                    n_panel_abs = n_panel_abs + 1;
                end
                ipc.te_last_indices(n_element) = n_panel_abs-1;
            end
            
            
            
            %% Join them together and formulate non-lifting problem
            
            % Now generate global panel record
            ipc.aps = airfoil_panel_set(panel_list);
            
            % Make source_induced_speeds_matrices
            [u_s_x_M, u_s_y_M] = source_induced_speeds_matrices(ipc.aps);
            % Convert to normal and tangential speeds
            [u_s_n_M] = normal_speeds_matrix(ipc.aps , u_s_x_M , u_s_y_M);
            
            % Make boundary condition normal velocities vector
            [ipc.u_n_BC_vector] = normal_BC_condition(ipc.aps);
            
            %% Vortex Part
            % Make source_induced_speeds_matrices
            [u_v_x_M, u_v_y_M] = vortex_induced_speeds_matrices(ipc.aps);
            % Convert to normal speeds
            [u_v_n_M] = normal_speeds_matrix(ipc.aps , u_v_x_M , u_v_y_M);
            
            % Make column vector for flow tangency condition (zero normal component)
            % Sum because all vortices have the same strenght
%             u_v_n_column = sum(u_v_n_M,2);
            
            % Now for multiple elements, recall that
%             sum(m,2) ==  sum(m(:,1:3),2) + sum(m(:,4:5),2)
            % So we can write
            % Allocate Space
            u_v_n_columns = zeros(size(u_v_n_M,1) , ipc.N_elements);
            for n_element = 1:ipc.N_elements
                range_tmp = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
                u_v_n_columns(:,n_element) = sum(u_v_n_M(:,range_tmp),2);
            end
            
            % Now join as usual
            ipc.u_n_M = [u_s_n_M , u_v_n_columns];
            A_tangency_lines = ipc.u_n_M;
            
            %% Now take care of kutta condition
            % Simple approach:
            %       impose equal pressure approximately
            %       by imposing same speed on the two panels adjacent to the trailing
            %       edge
            % Compute tangential components
            
            % This is the part that differs from the original formulation
            % for single element airfoils
            % Allocate Memory
            A_kutta_lines = zeros(ipc.N_elements , size(A_tangency_lines,2));
            for n_element = 1:ipc.N_elements
                [ipc.u_s_t_M] = tangential_speeds_matrix(ipc.aps , u_s_x_M , u_s_y_M);
                [ipc.u_v_t_M] = tangential_speeds_matrix(ipc.aps , u_v_x_M , u_v_y_M);
            
                % Extract element bounds for code readability
                f = ipc.te_first_indices(n_element);
                l = ipc.te_last_indices(n_element);                
                
                % Extract first and last lines for sources                
                u_s_t_line_first = ipc.u_s_t_M(f,:);
                u_s_t_line_last  = ipc.u_s_t_M(l,:);
            
                % Extract first and last lines for vortex components
                u_v_t_line_first = ipc.u_v_t_M(f,:);
                u_v_t_line_last  = ipc.u_v_t_M(l,:);
            
                % Now concatenate, (sum in 2nd direction, horizontally), as all vortices
                % have the same strength
                % Original
%                 u_v_t_concatenation_first = sum(u_v_t_line_first , 2);
%                 u_v_t_concatenation_last  = sum(u_v_t_line_last , 2);
                % Now multielement
                u_v_t_concatenation_first = zeros(size(u_v_t_line_first,1) , ipc.N_elements);
                u_v_t_concatenation_last  = zeros(size(u_v_t_line_last,1)  , ipc.N_elements);
                for n_element2 = 1:ipc.N_elements
                    range_tmp = ipc.te_first_indices(n_element2):ipc.te_last_indices(n_element2);
                    u_v_t_concatenation_first(:,n_element2) = sum(u_v_t_line_first(:,range_tmp),2);
                    u_v_t_concatenation_last(:,n_element2)  = sum(u_v_t_line_last(:,range_tmp) ,2);
                end
            
                % Now subtract tangential velocity coefficients of first and last panels to impose their
                % equality
                % Note: Actually we sum instead of subtracting because the tangential
                % versors of the upper and lower sides are oriented opposite to each other
                kutta_v_terms = u_v_t_concatenation_last + u_v_t_concatenation_first;
                kutta_s_terms = u_s_t_line_last + u_s_t_line_first;

                % And compose into matrix line
                A_kutta_lines(n_element,:) = [kutta_s_terms kutta_v_terms];
            end

            %% Now joint everything into big A matrix!
            ipc.A = [A_tangency_lines ; A_kutta_lines];
        end
        
        function generate_b_column(ipc)
            %% Now generate b column!
            
            % Start by make free stream components!
            % Make free stream speed vector (vertical)
            
            % Now account for free stream
            if ipc.rotation == false
                % For uniform free stream
                % Original (Winter 2010-2011)
                % ipc.u_inf_x_vec = - ipc.u_inf * ones(length(ipc.aps.panel_list) , 1) * cos(ipc.alpha);
                % ipc.u_inf_y_vec = - ipc.u_inf * ones(length(ipc.aps.panel_list) , 1) * sin(ipc.alpha);
                % Corrected (Spring 2017) (previous errors canceled out for most purposes)
                ipc.u_inf_x_vec =   ipc.u_inf * ones(length(ipc.aps.panel_list) , 1) * cos(ipc.alpha);
                ipc.u_inf_y_vec =   ipc.u_inf * ones(length(ipc.aps.panel_list) , 1) * sin(ipc.alpha);
            else
                % For rotating free stream
                [x_c_column, y_c_column] = centerpoint_locations_columns(ipc.aps);
                ipc.u_inf_x_vec = arrayfun(ipc.u_x_theta , x_c_column , y_c_column);
                ipc.u_inf_y_vec = arrayfun(ipc.u_y_theta , x_c_column , y_c_column);
            end
            
            % Rewrite it in normal and tangential modes
            [ipc.u_inf_n_vec] = normal_speeds_matrix(ipc.aps , ipc.u_inf_x_vec , ipc.u_inf_y_vec);
            [ipc.u_inf_t_vec] = tangential_speeds_matrix(ipc.aps , ipc.u_inf_x_vec , ipc.u_inf_y_vec);
            
            % First flow tangency condition
            b_tangency = (ipc.u_n_BC_vector - ipc.u_inf_n_vec);
            % Now kutta terms
            % % Original Single Element Version (Winter 2010-2011)
            % b_kutta = 0;
            % % Original multielement Version (verified) (Winter 2010-2011) (In use until 2017 update for this fork!)
            % b_kutta = zeros(ipc.N_elements , 1);
            % % New Single Element Version (Winter 2010-2011)
            % From conversation with Sachin (Spring 2015), I moved to this for b_kutta: U_inf<dot>tangent_1 - U_inf<dot>tangent_N 
            % Modification (re)introduced in Spring 2017 (not sure if I did it earlier for Jelmer, but can't find it back right now, right here!)
            % [u_inf_t_vec] = tangential_speeds_matrix(ipc.aps , ipc.u_inf_x_vec , ipc.u_inf_y_vec);
            % b_kutta = -( u_inf_t_vec(1)+u_inf_t_vec(end));
            % % New Multi Element Version (Spring 2017! Starts today!)
            %[ipc.u_inf_t_vec] = tangential_speeds_matrix(ipc.aps , ipc.u_inf_x_vec , ipc.u_inf_y_vec);
            b_kutta = -( ipc.u_inf_t_vec(ipc.te_first_indices)+ipc.u_inf_t_vec(ipc.te_last_indices));
                                    
            ipc.b = [ b_tangency ; b_kutta];            
        end      
        
        function solve(ipc)
            %% And Solve!
            ipc.q_vector = ipc.A \ ipc.b;
        end
        
        function post_process(ipc)
            % Post process lifting problem

            % Make tangential speeds matrix
            % Make tangential speed vortex contribution term
            % Original Version
%             u_v_t_column = sum(ipc.u_v_t_M,2);
            % Multi-element-version
            u_v_t_columns = zeros(size(ipc.u_v_t_M,1) , ipc.N_elements);
            for n_element = 1:ipc.N_elements
                element_indices = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
                u_v_t_columns(:,n_element) = sum(ipc.u_v_t_M(:,element_indices),2);
            end
            
            % Join into tangential coefficients matrix for sources and vortex strength
            ipc.u_t_M = [ipc.u_s_t_M , u_v_t_columns];
            
            % And get the speeds at the panels back
            ipc.u_t_vec = ipc.u_t_M * ipc.q_vector + ipc.u_inf_t_vec;
            
            % Compute Lift
            % From pressures
            % Define Cp
            ipc.cp = ipc.u_t_vec.^2;
            
            % Find normal versor directions
            [n_versor_x, n_versor_y] = normal_versors_column(ipc.aps);
            
            % Define force components on each panel
            f_x = ipc.cp .* n_versor_x;
            f_y = ipc.cp .* n_versor_y;
            
            % Now integrate using explicit trapeze rule
            % Get panel lengths
            l_column = panel_lengths_column(ipc.aps);
            
            % Obtain Net Forces
            % Original case
            ipc.F_x = sum(f_x .* l_column);
            ipc.F_y = sum(f_y .* l_column);
            % As there is no drag, all forces are lift            
            ipc.L_Cp = norm([ipc.F_x , ipc.F_y]);
            % Multielement case
            for n_element = 1:ipc.N_elements
                element_indices = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
                ipc.F_x(n_element) = sum(f_x(element_indices) .* l_column(element_indices));
                ipc.F_y(n_element) = sum(f_y(element_indices) .* l_column(element_indices));
            end
            ipc.L_Cp = -sin(ipc.alpha) * ipc.F_x + cos(ipc.alpha) * ipc.F_y;
            ipc.D_Cp =  cos(ipc.alpha) * ipc.F_x + sin(ipc.alpha) * ipc.F_y;
                                                
            % From circulation            
            for n_element = 1:ipc.N_elements
                element_indices = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
                ipc.L_gamma(n_element) = 2 * ipc.q_vector(end-ipc.N_elements + n_element) * sum(l_column(element_indices));
            end
            ipc.res2 = sqrt((ipc.L_Cp - ipc.L_gamma).^2);

            % And move on to obtain normal speeds
            u_n_vec = ipc.u_n_M * ipc.q_vector + ipc.u_inf_n_vec;
            ipc.res = sqrt(sum(u_n_vec.^2));
            
            ipc.cp_plot = ipc.u_t_vec.^2 - 1;            
        end
        
%% The Functions Below are only for postprocessing!        
        function [u_x_line, u_y_line] = induced_speed_line(ipc , x, y)
            % Get source and vortex influence coefficients at point
            [u_s_x_line, u_s_y_line] = ipc.aps.source_induced_speeds_line(x, y);
            [u_v_x_line, u_v_y_line] = ipc.aps.vortex_induced_speeds_line(x, y);
            % Now concatenate vortex line for multielement case
            u_v_x_line_concatenation = zeros(size(u_v_x_line,1) , ipc.N_elements);
            u_v_y_line_concatenation = zeros(size(u_v_y_line,1)  , ipc.N_elements);
            for n_element = 1:ipc.N_elements
                element_indices = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
                u_v_x_line_concatenation(:,n_element) = sum(u_v_x_line(:,element_indices),2);
                u_v_y_line_concatenation(:,n_element)  = sum(u_v_y_line(:,element_indices) ,2);
            end          
            % Join into complete influence line, compatible with solution
            % vector ipc.q_vector
            u_x_line = [u_s_x_line , u_v_x_line_concatenation];
            u_y_line = [u_s_y_line , u_v_y_line_concatenation];
        end
        
        function [u_x_M, u_y_M] = induced_speed_matrix(ipc , x, y)
            % Matrix for influence coefficients at list of points given in
            % x and y vectors
            
            N_points = length(x);
            % Allocate Memory
            u_x_M = zeros(N_points , length(ipc.q_vector));
            u_y_M = zeros(N_points , length(ipc.q_vector));
            parfor n_point = 1:N_points
                % Compute influence line for current point 
                [u_x_line, u_y_line] = induced_speed_line(ipc , x(n_point), y(n_point));
                % Stuff into matrices
                u_x_M(n_point , :) = u_x_line;
                u_y_M(n_point , :) = u_y_line;
            end
            
            % Done!
        end
        
%% These functions are for postprocessing and interaction with other solvers
        function [u, v] = induced_speed_on_many_points(ipc, x, y)
            % This is a grossly inneficient method, but it does the job
            % for now!
            
            % Serialize x and y input points (to deal with arbitrary dimension arrays)
            % x and y must have identical size!
            xV = x(:);
            yV = y(:);
            
            % Generate interaction lines
            [u_x_M, u_y_M] = induced_speed_matrix(ipc , xV, yV);
            
            % Compute Induced Speeds (to use all this, you still want to sum the free-stream!)
            uV = u_x_M * ipc.q_vector;
            vV = u_y_M * ipc.q_vector;
            
            % Now recast into suitable shape!
            % xM = magic(3); xV = xM(:); yV = 2*xV; yM = reshape(yV, size(xM)); yM/2 == xM; disp(['Success : ' num2str(isempty(find(not(yM/2 == xM))))]);
            u  = reshape(uV, size(x));
            v  = reshape(vV, size(x));  
        end
        
        function [x_c_column, y_c_column] = coupling_point_location_columns(ipc)
            % Should this really be the centerpoint? or the third quarter
            % chord point? to be checked before publishing! should do for
            % now (hopefully!) as long as we stay far enough for high
            % curvature regions! 
                [x_c_column, y_c_column] = centerpoint_locations_columns(ipc.aps);
        end
        
        function generate_solution_coupled(ipc, u_inf_x_vec_bis, u_inf_y_vec_bis)
            % Generates a solution for a coupled problem (called from dogoro/kirikou solvers)
            generate_b_column_coupled(ipc, u_inf_x_vec_bis, u_inf_y_vec_bis);
            ipc.solve();
        end
        
        function generate_b_column_coupled(ipc, u_inf_x_vec_bis, u_inf_y_vec_bis)
            %% Now generate b column for a coupled problem, called from
            
            % Here, external forcings (free-stream, other elements, etc, come as an input)
            
            % Rewrite forcings into normal and tangential components
            [ipc.u_inf_n_vec] = normal_speeds_matrix(ipc.aps , u_inf_x_vec_bis , u_inf_y_vec_bis);
            [ipc.u_inf_t_vec] = tangential_speeds_matrix(ipc.aps , u_inf_x_vec_bis , u_inf_y_vec_bis);
            
            % Write flow tangency condition
            b_tangency = (ipc.u_n_BC_vector - ipc.u_inf_n_vec);
            
            % Write kutta terms
            b_kutta = -( ipc.u_inf_t_vec(ipc.te_first_indices)+ipc.u_inf_t_vec(ipc.te_last_indices));
            
            % Compose b vector
            ipc.b = [ b_tangency ; b_kutta]; 
        end
            
    end
    
end

