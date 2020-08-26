classdef kirikou_single_actuator_solver < handle
    %KIRIKOU_SINGLE_ACTUATOR_SOLVER Manager Class in charge of initializing,
    %   running and postprocessing a single actuator (arbitrary loading)
    %   inviscid incompressible flow case including external vorticity
    %   supports
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                     %
    %   Kirikou -   A simple specialized 2d Vorticity Equation Solver for %
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
    
    %#ok<*NASGU>
    properties
        Ct                  =  8/9;     % Actuator Force Coefficient (C_F_a)

        x_v                 = -0.2;     % x-stance of lifting vortex pair (symmetry on x axis)
        r_v                 =  0.4;     % distance of lifting vortex pair to x axis (symmetry on x axis)
        gamma_v             =  0.4;     % strenght of each vortex in the pair

        % Free-Stream Charactheristics
        u_inf               = 1;                            % Unperturbed Wind Speed
        rho                 = 1;                            % Density

        % Rotor Inputs
        d1                  = 1;                            % Actuator Diameter
        f1_average                                          % Average force density over the actuator
        
        % X-Parametrization Inputs (single actuator streamtube)
        n_stances            =  20;                         % n_stances = n_panels + 1
        n_farwake_stances    =   0;                         % Add smooth far wake deformations
        n_filaments_per_side =   1;                         % Filaments per side

        % x-Parametrization Settings
        theta_start         = 0;                            % Zero means wake starts at actuator
        theta_end           = 0.9999*pi()/2;                % pi()/2 means wake extends down to infinity
        beta_vector_guess   = [1 1];                        % Coefficient vector (ones(...) means straight wake)
        xs                  = 0;                            % Origin of h-theta system in x-y plane
        ys                  = 0;                            %

        % y-Parametrization Settings
        lump_exp            = 1/5;                          % Lumping coefficient for y-discretization
        
        % Solver Parameters (Stopping Criteria)
        information_period  = 400;                          % Display information/plot every information_period iterations
        max_RES_stretch     = 1e-12;                        % Maximum Stretching (extension) residual (vorticity conservation)
        max_RES_shape       = 1e-12;                        % Maximum Shape residual (volume flow accross segments)
        % Solver Parameters (Relaxation) (0.005 0.005 0.0)
        relax = 0.005;                                      % Stretching Relaxation
        dt    = 0.005;                                      % Shape (Virtual) Time Step
        dtx   = 0.0;                                        % Shape (Virtual) Time Step (x movement)
        % Computational Parameters
        use_parallel        = false;                        % If true, compute induction in parallel        
        N_iterations = 40000;                               % Total number of iterations for solving wake
        
        % % Processed Parameters
        r1                                                  % Radius of  Upstream Actuator
        S1                                                  % Surface of Upstream Actuator
        f_w                                                 % Loading function
        y_release                                           % Filament Release Points
        delta_f_w                                           % Loading Steps that Generate Filaments (like a discrete derivative)
        n_panels                                            % Number of Panels on Mainwake
        n_filaments                                         % Total number of filaments

        % % Solution Parameters
        % Wake Geometry
        x_start                                             % List of Panel Startpoints 
        y_start                                             % 
        x_end                                               % List of Panel Endpoints
        y_end                                               %
        % Wake Strenght
        gamma0                                              % Reference Circulation Densities
        gamma                                               % Actual    Circulation Densities
        % Speeds (Induced+FreeStream) on Wake
        s_l                                                 % Tangential Speed Component on Panels
        t_n                                                 % Normal     Speed Component on Panels
        
        % Index Arrays (used for boundary condition and fudging enforcement on generic rotors, including inhomogeneous ones)
        id_all                                              % List all panel edge points
        id_no_end                                           % List all panel edge points but the ones at filament end (very far dowsntream)
        id_end                                              % List filament end points (very far dowsntream) indices in complete array of panel edge points 
        id_start                                            % List filament start points (on actuator plane) indices in complete array of panel edge points 
        
        % Global Residual Values
        RES_stretch                                         % Wake vorticity strenght residual 
        RES_shape                                           % Wake shape misalignment residual
        
        % Residual convergence logging arrays
        log_n_array                                         % Log storage iterations
        log_RES_shape_array                                 % Log shape residual
        log_RES_stretch_array                               % Log stretch (vorticity) residual
        
        % % Postprocessing
        % Machine Parameters
        x1                                                  % x Position        (First  Actuator )
        a1                                                  % Induction Factor  (Absolute)
        m1                                                  % Mass Flux         (First  Actuator )
        p1                                                  % Power             (First  Actuator )
        Cp1                                                 % Power Coefficient (First  Actuator )
        % Complete Machine
        p                                                   % Power of Complete Machine
        Cp                                                  % Power Coefficient of Complete Machine
        
        % Velocity Fields
        x_mesh                                              % Mesh of x positions to compute velocity field (for plotting only)
        y_mesh                                              % Mesh of y positions to compute velocity field (for plotting only)
        u_mesh                                              % Velocity field (x-component) on mesh          (for plotting only)
        v_mesh                                              % Velocity field (y-component) on mesh          (for plotting only)
        u_norm_mesh                                         % Velocity field (magnitude)   on mesh          (for plotting only)
        
        % % Handles to Other Objects
        VS                                                  % Handle to vortex element object
        DL                                                  % Handle to disk loading object
        
        % % Results Bundle
        RES                                                 % Bundles key machine parameters for compact storage to explore performance trends

    end
    
    methods
        function DSD = kirikou_single_actuator_solver(Ct, x_v, r_v, gamma_v)
            % Constructor Method for the Single Actuator Solver Class
            % % Import Actuator Parameters
            % Actuator Thrust Coefficient (imposes force density)
            DSD.Ct = Ct;
            % Vortex pair streamwise stance
            DSD.x_v = x_v;
            % Vortex pair height
            DSD.r_v = r_v;
            % Circulation of each vortex in pair
            DSD.gamma_v = gamma_v;
        end
        
        function preprocess_run_postprocess(DSD)
            % Prepocess Inputs
            DSD.preprocess_inputs();
            % Discretize with Straight Wake Tube as Initial Guess
            DSD.define_discretization_and_initial_geometry();
            % Define Reference Panel Strenghts
            DSD.define_reference_panel_strenghts();
            % Now formulate Initial Guesses
            DSD.formulate_initial_guesses();
            % Create Vortex Element Object (with state provided by initial guesses)
            DSD.create_vortex_element_object();
            % Now Start Solver
            DSD.run_solver();
            % Now compute Machine Parameters
            DSD.compute_machine_parameters();
            % Generate Velocity Fields
            DSD.generate_velocity_fields();
        end
        
        function preprocess_inputs(DSD)
            % Compute Actuator Radius
            DSD.r1 = DSD.d1/2;
            % Set surface (2d surface is just diameter/width)
            DSD.S1 = DSD.d1;
            
            % Set force from Ct
            F1         = - 0.5*DSD.rho*DSD.S1*(DSD.u_inf^2)*DSD.Ct;   % Actuator Force 
            DSD.f1_average =   F1 / DSD.d1;                  % Average Actuator Force per unit radius
            % From target terminal wake velocity
            % ue_target    = 1/3
            % f1_average   = - 0.5 * rho * u_inf.^2 * (1-ue_target^2);
            
            %% Define Force Function (New model)
            % % New Model
            DSD.DL = disk_loading(DSD.f1_average);
            % Loading Function
            DSD.f_w = @(y_over_r) DSD.DL.f_w_function(y_over_r);
            % Original (did not fall to zero on sides
            %f_w = @(y_over_r) C0(y_over_r) .* f_w_disk_flag(y_over_r) * f1_average;
        end
        
        function define_discretization_and_initial_geometry(DSD)
            %% Geometry definition (x-Parametrization)
            % % Generate geometry and stance steps
            % Addapt
            beta_vector = DSD.r1*DSD.beta_vector_guess;                          % Scale Coefficients to radius
            
            % % Generate range of theta values for panel definition
            % (almost always between 0 and pi()/2)
            theta_range = linspace(DSD.theta_start, DSD.theta_end, DSD.n_stances-DSD.n_farwake_stances);  
            % Recast theta range as column vector (so pas tosses)
            theta_range = theta_range(:);                               
            
            % Generate h values from theta values (beta_vector sets initial
            % wake shape guess, currently almost unused, best results are
            % found with straight line as initial guess, this code stays
            % here to allow later improvements) 
            [h , ~] = parametrization_functions.h_streamline(theta_range , beta_vector);
            
            % Map h-theta array to x-y plane (map [0,pi/2] interval
            % into [0,+inf] interval
            [x,y, ~, ~, ~, ~] = parametrization_functions.S_transformation(h, theta_range, DSD.xs, DSD.ys);
            disp(['Wake Lenght : ' num2str(max(x)/DSD.d1) ' diameters']);
            
            % % Add log spaced stances for better far wake resolution
            %   Improves mass conservation, but makes convergence more
            %   difficult!
            if DSD.n_farwake_stances > 0
                % Remove last stance (usually too far)
                x = x(1:end-1);
                x_far_wake = x(end) * logspace(0, 2, DSD.n_farwake_stances+2);
                x = [x(1:end-1) ;  x_far_wake(:)];
                y = DSD.r1 * ones(size(x));
            end
                                    
            %% Geometry definition (Construct generating filament)
            % % Generate discretization from initial geometry/stance steps
            % Determine number of panels
            DSD.n_panels = DSD.n_stances - 1;
            
            % Now segment stances into sets of reference panels
            x_start_generating = x(1:(end-1));
            y_start_generating = y(1:(end-1));

            x_end_generating   = x(2:end);
            y_end_generating   = y(2:end);

            % Recast generating filament vectors into columns
            x_start_generating = x_start_generating(:);
            y_start_generating = y_start_generating(:);

            x_end_generating   = x_end_generating(:);
            y_end_generating   = y_end_generating(:);
            
            %% Geometry definition (y-Discretization)
            % Segment y Direction in Constant Loading Steps
            y_over_r_side_outer = linspace(1/DSD.n_filaments_per_side,  1, DSD.n_filaments_per_side).^(DSD.lump_exp);
            y_over_r_side_inner = [0 , y_over_r_side_outer(1:(DSD.n_filaments_per_side-1))];
            
            % Determine force step for each filament (proportional to vorticity generation rate)
            delta_f_w_side = DSD.f_w(y_over_r_side_outer) - DSD.f_w(y_over_r_side_inner);
            
            %% Geometry definition (Replicate Filaments over Upper Half Side)
            % Allocate Panel Position Vectors
            x_start_upper_side   = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            y_start_upper_side   = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);

            x_end_upper_side     = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            y_end_upper_side     = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
           
            % Allocate Unperturbed Strenght Vectors
            gamma0_upper_side     = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            
            % Allocate Filament Index Vectors (id_(...))
            id_start_upper_side  = zeros(DSD.n_filaments_per_side, 1);
            id_end_upper_side    = zeros(DSD.n_filaments_per_side, 1);
            
            % Allocate y Release reference array
            y_release_upper_side = zeros(DSD.n_filaments_per_side, 1);
            delta_f_w_upper_side = zeros(DSD.n_filaments_per_side, 1); 
            
            % Fill side arrays using generating filament arrays ((...)_generating)
            for n_filament = 1:DSD.n_filaments_per_side
                % Compute release position of filament (reverse to move from edge to
                % center)
                y_release_upper_side(n_filament) = y_over_r_side_outer(DSD.n_filaments_per_side - n_filament +1) * DSD.r1;
                % Compute force step generating the filament (reverse to move from edge
                % to center)
                delta_f_w_upper_side(n_filament) = delta_f_w_side(DSD.n_filaments_per_side - n_filament +1);
                
                % Set indexes
                id_start_upper_side(n_filament) = 1 + (n_filament-1)*DSD.n_panels;
                id_end_upper_side(  n_filament) =     (n_filament  )*DSD.n_panels;
                
                % Fill side panel arrays using generating filament and (just/last
                % generated) index entries
                x_start_upper_side(id_start_upper_side(n_filament):id_end_upper_side(n_filament)) = ...
                    x_start_generating;
                y_start_upper_side(id_start_upper_side(n_filament):id_end_upper_side(n_filament)) = ...
                    y_start_generating * (y_release_upper_side(n_filament) / DSD.r1);
                x_end_upper_side(  id_start_upper_side(n_filament):id_end_upper_side(n_filament)) = ...
                    x_end_generating;
                y_end_upper_side(  id_start_upper_side(n_filament):id_end_upper_side(n_filament)) = ...
                    y_end_generating   * (y_release_upper_side(n_filament) / DSD.r1);
                
                % Fill Unperturbed Strenght Vector (needs modification for yawed flows!)
                gamma0_upper_side( id_start_upper_side(n_filament):id_end_upper_side(n_filament)) = ...
                    delta_f_w_upper_side(n_filament) /  DSD.u_inf;
            end
            
            %% Geometry definition (Replicate Filaments over Lower Half Side)
            % Allocate Panel Position Vectors
            x_start_lower_side   = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            y_start_lower_side   = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            
            x_end_lower_side     = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            y_end_lower_side     = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            
            % Allocate Unperturbed Strenght Vectors
            gamma0_lower_side     = zeros(DSD.n_panels*DSD.n_filaments_per_side, 1);
            
            % Allocate Filament Index Vectors (id_(...))
            id_start_lower_side  = zeros(DSD.n_filaments_per_side, 1);
            id_end_lower_side    = zeros(DSD.n_filaments_per_side, 1);
            
            % Allocate y Release reference array
            y_release_lower_side = zeros(DSD.n_filaments_per_side, 1);
            delta_f_w_lower_side = zeros(DSD.n_filaments_per_side, 1);
            
            % Fill side arrays using generating filament arrays ((...)_generating)
            for n_filament = 1:DSD.n_filaments_per_side
                % Compute release position of filament (direct to move from center to
                % edge)
                y_release_lower_side(n_filament) = - y_over_r_side_outer(n_filament) * DSD.r1;
                % Compute force step generating the filament (direct to move from
                % center to edge)
                delta_f_w_lower_side(n_filament) =   delta_f_w_side(n_filament);
                
                % Set indexes
                id_start_lower_side(n_filament) = 1 + (n_filament-1)*DSD.n_panels;
                id_end_lower_side(  n_filament) =     (n_filament  )*DSD.n_panels;
                
                % Fill side panel arrays using generating filament and (just/last
                % generated) index entries
                x_start_lower_side(id_start_lower_side(n_filament):id_end_lower_side(n_filament)) = ...
                    x_start_generating;
                y_start_lower_side(id_start_lower_side(n_filament):id_end_lower_side(n_filament)) = ...
                    y_start_generating * (y_release_lower_side(n_filament) / DSD.r1);
                x_end_lower_side(  id_start_lower_side(n_filament):id_end_lower_side(n_filament)) = ...
                    x_end_generating;
                y_end_lower_side(  id_start_lower_side(n_filament):id_end_lower_side(n_filament)) = ...
                    y_end_generating   * (y_release_lower_side(n_filament) / DSD.r1);
                
                % Fill Unperturbed Strenght Vector (inverted for symmetry!)(needs modification for yawed flows!)
                gamma0_lower_side( id_start_lower_side(n_filament):id_end_lower_side(n_filament)) = ...
                    - delta_f_w_lower_side(n_filament) /  DSD.u_inf;
            end
            
            %% Geometry definition (Join Sides)
            % Set Total number of Filaments
            DSD.n_filaments = DSD.n_filaments_per_side * 2;
            
            % Join Panel Position Vectors
            DSD.x_start = [x_start_upper_side(:) ; x_start_lower_side(:)];
            DSD.y_start = [y_start_upper_side(:) ; y_start_lower_side(:)];
            
            DSD.x_end   = [  x_end_upper_side(:) ; x_end_lower_side(:)  ];
            DSD.y_end   = [  y_end_upper_side(:) ; y_end_lower_side(:)  ];
            % Join Unperturbed Str%HELPER_FUNCTIONS is a simple class collecting static helper functions
    %for the Dogoro and Kirikou solversenght Vectors
            DSD.gamma0  = [ gamma0_upper_side(:) ; gamma0_lower_side(:) ];
            
            % Now Allocate and Join Filament Index Vectors (id_(...))
            % For Filament Starts
            DSD.id_start = zeros(DSD.n_filaments, 1);
            DSD.id_start(1:DSD.n_filaments_per_side)     = id_start_upper_side;
            DSD.id_start(DSD.n_filaments_per_side+1:end) = id_start_lower_side + DSD.n_panels * DSD.n_filaments_per_side;
            
            % For Filament Ends
            DSD.id_end = zeros(DSD.n_filaments, 1);
            DSD.id_end(1:DSD.n_filaments_per_side)     = id_end_upper_side;
            DSD.id_end(DSD.n_filaments_per_side+1:end) = id_end_lower_side     + DSD.n_panels * DSD.n_filaments_per_side;
            
            % And finally fill y Release reference arrays
            DSD.y_release = [y_release_upper_side(:) ; y_release_lower_side(:)];
            DSD.delta_f_w = [delta_f_w_upper_side(:) ; delta_f_w_lower_side(:)];
            
            
        end
        
        function define_reference_panel_strenghts(DSD) %#ok<MANU>
%             % % % Initial guess                        
%             
%             % % Generate array of forces
%             f_panels = ones(1, DSD.n_panels);
%             % Make early wake, due to 1st actuator only
%             f_panels(1:DSD.n_nearwake_stances) = DSD.f1;
%             % Make late wake, superposition of 1st and 2nd actuator
%             f_panels(DSD.n_nearwake_stances+1:DSD.n_panels) = DSD.f1 + DSD.f2;
%             
%             % % Generate array of reference circulation densities
%             % The reference circulation density of the panels is based on
%             % the situation in which the tangential speed on the panel
%             % equals free-stream.
%             %
%             % The equal circulation density is found by the solver, which
%             % scaling this reference until the extension (stretch) residual
%             % is minimized 
%             % 
%             % Upper side:
%             DSD.gamma0(1:DSD.n_panels)       = f_panels./DSD.u_inf;
%             % Lower side:                         (symmetric of upper side)
%             DSD.gamma0((DSD.n_panels+1):(2*DSD.n_panels)) = ...
%                                             - DSD.gamma0(1:DSD.n_panels);
%             % Turn into Column Vector
%             DSD.gamma0 = DSD.gamma0(:);
                        
        end
        
        function formulate_initial_guesses(DSD)
            %% Complete Memory Allocation and Initialization
            % Allocate arrays for tangential (s_l) and normal (t_n) speeds components
            % s_l = zeros(n_panels*n_filaments, 1);
            DSD.t_n = zeros(DSD.n_panels*DSD.n_filaments, 1);
            
            % Formulate initial guess for vortex strenght
            DSD.gamma = DSD.gamma0;
            
            % Make an Index Vector to index all panels
            DSD.id_all = 1:(DSD.n_panels*DSD.n_filaments);
            
            % Make an Index Vector to index all panels except the last ones of each filament
            % Initialize with full panel list
            DSD.id_no_end = 1:(DSD.n_panels*DSD.n_filaments);
            % Remove elements
            for n_filament = 1:DSD.n_filaments
                DSD.id_no_end(DSD.id_no_end == DSD.id_end(n_filament)) = [];
            end
            % Axi: Remove second round of elements
            % Remove elements
            for n_filament = 1:DSD.n_filaments
                DSD.id_no_end(DSD.id_no_end == (DSD.id_end(n_filament)-1)) = [];
            end
        end
        
        function create_vortex_element_object(DSD)
            % Create vortex segment object
            if not(DSD.use_parallel)
                DSD.VS = constant_strenght_vortex_segment_2d(DSD.x_start, DSD.y_start, DSD.x_end, DSD.y_end, DSD.gamma);
            else
                DSD.VS = constant_strenght_vortex_segment_2d_parallel(DSD.x_start, DSD.y_start, DSD.x_end, DSD.y_end, DSD.gamma);
            end
 
            % Create vortex pair object
            VP = singular_vortex_pair_2d(DSD.x_v, DSD.r_v, DSD.gamma_v);
            % Hook to vortex element object
            DSD.VS.VP = VP;
        end
        
        function run_solver(DSD)
            figure(9); hold on; grid on;
            
            % Allocate Residual convergence logging arrays
            DSD.log_n_array           = zeros(1, DSD.N_iterations);
            DSD.log_RES_shape_array   = zeros(1, DSD.N_iterations);
            DSD.log_RES_stretch_array = zeros(1, DSD.N_iterations);
            
            for n=1:DSD.N_iterations
                % % Compute Movement on Panel Centers
                %   Zero-out x movement (DSD.dtx=0) makes convergence a lot
                %   slower but ensures that we preserve a very fine
                %   discretization in critical areas.
                dx_center = DSD.t_n.* DSD.VS.x_n_unit_vector * DSD.dtx;
                dy_center = DSD.t_n.* DSD.VS.y_n_unit_vector * DSD.dt;
                
                % New: bypass deformation of last (most downstream) panel to prevent
                % wake from closing on itself due to starting vortex, this is an early
                % approach, can be refined with semi-infinite vortex segment later on!
                dx_center(DSD.id_end) = 0;
                dy_center(DSD.id_end) = 0;
                
                % Axi: bypass deformate on two last (most downstream)
                % panels
                dx_center(DSD.id_end-1) = 0;
                dy_center(DSD.id_end-1) = 0;
                
                % Filament by Filament Operations (Move Panels and Update Connectivity)
                for n_filament = 1:DSD.n_filaments
                    % % Move Panels
                    DSD.x_end(DSD.id_start(n_filament):DSD.id_end(n_filament)) = ...
                        DSD.x_end(DSD.id_start(n_filament):DSD.id_end(n_filament)) + cumsum(dx_center(DSD.id_start(n_filament):DSD.id_end(n_filament)));
                    DSD.y_end(DSD.id_start(n_filament):DSD.id_end(n_filament)) = ...
                        DSD.y_end(DSD.id_start(n_filament):DSD.id_end(n_filament)) + cumsum(dy_center(DSD.id_start(n_filament):DSD.id_end(n_filament)));
                    % % Synchronize starts with ends
                    DSD.x_start(DSD.id_start(n_filament)+1:DSD.id_end(n_filament)) =  DSD.x_end(DSD.id_start(n_filament):DSD.id_end(n_filament)-1);
                    DSD.y_start(DSD.id_start(n_filament)+1:DSD.id_end(n_filament)) =  DSD.y_end(DSD.id_start(n_filament):DSD.id_end(n_filament)-1);
                end
                
                % Update Vortex System Object with New Wake Shape and
                % pannel strenghts
                DSD.VS.update_independent_fields(DSD.x_start, DSD.y_start, DSD.x_end, DSD.y_end, DSD.gamma);
               
                % Compute Induced Speeds
                %   Speeds induced by wake, on panel centerpoints
                [~ , u_ctr , v_ctr ] = induced_speed_on_many_points(DSD.VS, DSD.VS.x_center, DSD.VS.y_center);
                % Obtain Tangential Speed Component and add Free-Stream
                DSD.s_l = DSD.u_inf .* DSD.VS.x_l_unit_vector + (u_ctr .* DSD.VS.x_l_unit_vector + v_ctr .* DSD.VS.y_l_unit_vector);
                % Obtain Normal Speed Component and add Free-Stream
                DSD.t_n = DSD.u_inf .* DSD.VS.x_n_unit_vector + (u_ctr .* DSD.VS.x_n_unit_vector + v_ctr .* DSD.VS.y_n_unit_vector);
                
                % Filter Speeds to Compute Residuals (Exclude last panels!)
                t_n_RES     = DSD.t_n( DSD.id_no_end);
                l_RES       = DSD.VS.l(DSD.id_no_end);
                
                % Compute Global Residuals (Original)
                %DSD.RES_stretch = sum((DSD.gamma - DSD.gamma0./DSD.s_l).^2 .* DSD.VS.l);
                %DSD.RES_shape   = sum(t_n_RES.^2 .* l_RES);
                % Compute Global Residuals (Axi)
                DSD.RES_stretch = sum((DSD.gamma(DSD.id_no_end) - DSD.gamma0(DSD.id_no_end)./DSD.s_l(DSD.id_no_end)).^2 .* DSD.VS.l(DSD.id_no_end));
                DSD.RES_shape   = sum(t_n_RES.^2 .* l_RES);
                
                
                % Find local stretch residual, to obtain  RES_stretch derivatives to
                % gamma vector components
                % 2d, Planar flow
                RESl_stretch        =     (DSD.gamma - DSD.gamma0 ./ DSD.s_l                        ).^2 .* DSD.VS.l;
                RESl_stretch_dgamma = 2 * (DSD.gamma - DSD.gamma0 ./ DSD.s_l                        )    .* DSD.VS.l;
                % Non adaptive (axi-symmetric ?)
                %DSD.RES_stretch     = 0;
                %RESl_stretch        = 0;
                %RESl_stretch_dgamma = 1;
                
                % Axisymmetric flow
                %RESl_stretch        =     (DSD.gamma - DSD.gamma0 .* sqrt(abs(DSD.s_l(1)./DSD.s_l)) ).^2 .* DSD.VS.l;
                %RESl_stretch_dgamma = 2 * (DSD.gamma - DSD.gamma0 .* sqrt(abs(DSD.s_l(1)./DSD.s_l)) )    .* DSD.VS.l;
                
                % Force fixed gamma near end (last and penultimate panels)
                DSD.gamma(DSD.id_end)   = DSD.gamma(DSD.id_end-2);
                DSD.gamma(DSD.id_end-1) = DSD.gamma(DSD.id_end-2);
                
                % Update gamma
                DSD.gamma               = DSD.gamma - DSD.relax * RESl_stretch ./ RESl_stretch_dgamma;
                % Keep filaments without circulation at 0 (correct NaN appearance due to zero RESl_stretch_dgamma for 0 gamma!)
                DSD.gamma(DSD.gamma0==0) = 0;
                
                % Display iteration information if demanded
                if or(n==1, round(n/DSD.information_period) == n/DSD.information_period)
                    DSD.display_periodic_information(n);
                end
                
                % Save Residuals to Loop
                DSD.log_n_array(n)           = n;
                DSD.log_RES_shape_array(n)   = DSD.RES_shape;
                DSD.log_RES_stretch_array(n) = DSD.RES_stretch;
                
                % Stop Loop if Residuals below ultimate tolerance
                if and(DSD.RES_stretch < DSD.max_RES_stretch , DSD.RES_shape < DSD.max_RES_shape)
                    disp('Solution Complete: Residual targets reached!')
                    DSD.display_periodic_information(n);                    
                    break
                end
                
            end
        end
        
        function display_periodic_information(DSD, n)
            % Write some stuff out!
            disp(['iter ', num2str(n) , '   RES_stretch = ' , num2str(DSD.RES_stretch ,'%8.4e') , ...
                '   RES_shape = ' , num2str(DSD.RES_shape , '%8.4e')        ]);
                %'   RES_strecht_dgamma = ' , num2str(sum(RESl_stretch_dgamma) , '%8.4e')]);
                
            % Plot Evolution of Residuals!
            figure(9)
            plot(n, log10(DSD.RES_shape), 'xk', n, log10(DSD.RES_stretch), 'ok');
            xlabel('N_{iter}'); ylabel('Residuals'); grid on; hold on;
            legend('Res_{shape}', 'Res_{\gamma}' );
            
            figure(11)
            plot(DSD.x_start, DSD.y_start, 'o-'); hold on
            plot(DSD.x_start, DSD.gamma, 'o-'); hold off
            xlabel('x (z)'); ylabel('y (r)'); grid on;
            
            drawnow;
            disp('New iteration')
        end
        
        function compute_machine_parameters(DSD)            
            %% Postprocessing (Reconstruct force distribution)
            % Compute force on patches from force steps!
            f1_actual_patches = - cumsum(DSD.delta_f_w .* sign(DSD.y_release));
            % Reinterpolate through it with nearest neighbour interpolar
            f1_actual = @(y) interp1(DSD.y_release, f1_actual_patches , y, 'next', 0);
            
            %% Postprocessing (Compute Induction and Power Coefficient of Actuator) (New Method)
            DSD.x1       = 0;
            eps_offset   = 1e-6;
            y_range1 = linspace(-(DSD.r1-eps_offset), (DSD.r1-eps_offset), 10000);
            x_range1 = DSD.x1 * ones(size(y_range1));
            % Compute wake effects on main extractor (all wakes!)
            [~ , u_actuator1   , ~]   = induced_speed_on_many_points(DSD.VS, x_range1 , y_range1);
            % Compute normal speed on actuator (only x component is needed,
            % as we would take dot product of complete vector with actuator normal
            % which is aligned with x axis!)
            u_a1 = u_actuator1 + DSD.u_inf;
            % Now compute induction factor
            % a_actuator   = mean(u_inf-u_a_actuator(not(isnan(u_a_actuator))));
            DSD.a1   = trapz(y_range1 , DSD.u_inf - u_a1) ./ (DSD.d1 *DSD.u_inf);
            % Compute mass flow on surface
            %m1   = trapz(y_range1, u_a1);
            % And perform integral for power!
            DSD.p1   = trapz(y_range1, f1_actual(y_range1) .* u_a1);
            % Reach power coefficient!
            DSD.Cp1  = DSD.p1 / (0.5*DSD.rho*DSD.S1*(DSD.u_inf^3));                   % Finally reach Cp
            
            Cp_num = DSD.Cp1;
            ua_num = trapz(y_range1, u_a1);
            a1_num = DSD.a1;
            
            %% Compare with Theorethical Values
            % Compute induced speeds on vortex pair (upper side)
            [~, u_v, v_v] = DSD.VS.induced_speed_on_many_points(DSD.VS.VP.x_v, DSD.VS.VP.r_v);
            u_v = u_v + 1;
            % Find norm of velocity on upper vortex of pair
            u_mag_v = sqrt(u_v^2 + v_v^2);
            
            % Find lift of upper vortex from pair
            L_v = DSD.rho*u_mag_v*DSD.VS.VP.gamma_v;
            % Find x component of force on upper vortex
            F_v_x = v_v / u_mag_v * L_v;
            %F_v_y = u_v / u_mag_v * L_v;
            
            % Double, because there are two vortices
            F_b = F_v_x*2;
            % Transform into coefficient
            C_F_b = F_b / (1/2 * DSD.rho * DSD.u_inf^2 * DSD.S1);
            
            % Theoretical Functions
            ue_fun      = @(C_f_a) sqrt(1 - C_f_a);
            % Cp_fun      = @(ue, xi) (1/2)*xi.*(ue + 1) + (1/2)*(ue + 1).*(ue.^2 - 1);
            Cp_fun      = @(ue, C_F_b) (1/2).*(ue + 1).*((ue.^2 - 1) - C_F_b);
            
            % Theoretical Values
            C_F_a       = DSD.Ct;                           % Ct as defined in older notation, for actuator forces only
            %ue_theo     = ue_fun(C_F_a);
            Cp_theo     = Cp_fun(ue_fun(C_F_a), C_F_b);
            ua_theo     = 1/2 * (1 + ue_fun(C_F_a)) - 1/2* C_F_b / (ue_fun(C_F_a) - 1);
            a1_theo     = 1 - ua_theo/DSD.u_inf;
            
            %% Display Diagnostic Text
            % Inputs
            disp('-- Problem Inputs --')
            disp(['r_a      = ' , num2str(DSD.r1) ]);
            disp(['phi_a    = ' , num2str(DSD.f1_average) ]);
            disp(['x_v      = ' , num2str(DSD.VS.VP.x_v) ]);
            disp(['r_v      = ' , num2str(DSD.VS.VP.r_v) ]);
            disp(['gamma_v  = ' , num2str(DSD.VS.VP.gamma_v) ]);
            % Force
            disp('-- Force Coefficient --')
            disp(['C_F_a    = ' , num2str(C_F_a) ]);
            disp(['C_F_b    = ' , num2str(C_F_b) ]);
            % Power Coefficient
            e_Cp = (Cp_num - Cp_theo) / Cp_theo;
            disp('-- Power Coefficient --')
            disp(['Cp_num   = ' , num2str(Cp_num) ]);
            disp(['Cp_theo  = ' , num2str(Cp_theo)]);
            disp(['e_Cp     =  ' , num2str(e_Cp*100)  , '%']);
            % Speed on Actuator
            e_ua = (ua_num - ua_theo) / ua_theo;
            disp('-- Average Normal Speed on Actuator --')
            disp(['ua_num   = ' , num2str(ua_num) ]);
            disp(['ua_theo  = ' , num2str(ua_theo)]);
            disp(['e_ua     = ' , num2str(e_ua*100)  , '%']);
            % Induction Coefficients
            e_a1 = (a1_num - a1_theo) / a1_theo;
            disp('-- Induction Factor --')
            disp(['a1_num   = ' , num2str(a1_num) ]);
            disp(['a1_theo  = ' , num2str(a1_theo)]);
            disp(['e_a1     = ' , num2str(e_a1*100)  , '%']);
            % Solution Parameters
            disp('-- Solution Residuals --')
            disp(['Strecht  = ' , num2str(DSD.RES_stretch) ]);
            disp(['Shape    = ' , num2str(DSD.RES_shape)   ]);
            
            %% Save to Results Bundle
            DSD.RES.r_a         = DSD.r1;
            DSD.RES.phi_a       = DSD.f1_average;
            DSD.RES.x_v         = DSD.x_v;
            DSD.RES.r_v         = DSD.r_v;
            DSD.RES.gamma_v     = DSD.gamma_v;
            
            DSD.RES.C_F_a       = C_F_a;
            DSD.RES.C_F_b       = C_F_b;
            
            DSD.RES.Cp_num      = Cp_num;
            DSD.RES.Cp_theo     = Cp_theo;
            DSD.RES.e_Cp        = e_Cp;
            
            DSD.RES.ua_num      = ua_num;
            DSD.RES.ua_theo     = ua_theo;
            
            DSD.RES.a1_num      = a1_num;
            DSD.RES.a1_theo     = a1_theo;
            
            DSD.RES.RES_stretch = DSD.RES_stretch;
            DSD.RES.RES_shape   = DSD.RES_shape;
            DSD.RES.VS          = DSD.VS;
            
            %%
%           % Original code from dogoro double solver, more generic,
%             move to this approach when time allows
%             % Find Position and Radius of First Actuator
%             DSD.x1 = DSD.VS.x_start(1);
%             DSD.r1 = DSD.VS.y_start(1);
%             DSD.d1 = DSD.r1*2;
%             DSD.S1 = DSD.d1;
%             % Compute induction, mass flux and power on surface of second actuator
%             [DSD.a1 , DSD.m1, DSD.p1] = DSD.compute_induction_on_cross_section(DSD.x1, DSD.r1, DSD.f1);
%             % And finally reach power coefficient!
%             DSD.Cp1  = DSD.p1 / (0.5*DSD.rho*DSD.S2*(DSD.u_inf^3));


             
        end
        
        function [a, m, p] = compute_induction_on_cross_section(DSD, x_section, r_section, f_section)
            % Find Section Diameter and Surface
            d_section = 2*r_section;            
            % Discretize the Cross-Section
            y_range = linspace(-(r_section-eps), (r_section-eps), 1000);
            x_range = x_section * ones(size(y_range));
            % Compute wake effects on main extractor (all wakes!)
            [~ , u_section   , ~]   = DSD.VS.induced_speed_on_many_points(x_range , y_range);
            % Compute normal speed on section (only x component is needed,
            % as we would take dot product of complete vector with actuator normal
            % which is aligned with x axis!)
            u_a = u_section + DSD.u_inf;
            % Now compute induction factor  (trapezoid rule)
            a   = trapz(y_range , DSD.u_inf - u_a) ./ (d_section *DSD.u_inf);
            % Compute mass flow on surface  (trapezoid rule)
            m   = trapz(y_range, u_a);
            % Integrate power!              (trapezoid rule)
            p   = trapz(y_range, f_section .* u_a);
        end
        
        function generate_velocity_fields(DSD)
            
            % % Define Velocity Field Generation Window
            % Bounds
            x_min = -4; x_max = 4;
            y_min = -2; y_max = 2;
            % Resolution
            x_range = linspace(x_min, x_max, 1000);
            y_range = linspace(y_min, y_max, 400);
            
            % % Make Mesh For Point Computation
            [DSD.x_mesh,DSD.y_mesh] = meshgrid(x_range, y_range);
            
            % % Compute Induced Speeds on Mesh Points 
            % (there is an error/instability on the potential computation,
            % but the speed components are correct)
            [~ , u_induced_mesh , v_induced_mesh] = induced_speed_on_many_points(DSD.VS, DSD.x_mesh, DSD.y_mesh);
            
            % Add Free-Stream to Induced Speeds
            DSD.u_mesh = u_induced_mesh + DSD.u_inf;
            DSD.v_mesh = v_induced_mesh;
            % And compute velocity norm
            DSD.u_norm_mesh = sqrt(DSD.u_mesh.^2 + DSD.v_mesh.^2);
        end
        
        function plot_velocity_field_with_streamlines(DSD, varargin)
            %% Plot Speed Magnitude with Streamlines
            % Get to right figure
            if isempty(varargin)
                figure(4); hold on
            else
                try 
                    axes(varargin{1})
                catch
                    try
                        figure(varargin{1});
                    catch
                    end
                end
            end
            
            % Plot actuation lines
            plot3([DSD.x1 DSD.x1], [-DSD.r1 DSD.r1] , 3*[1 1], '-' , 'Color', [0 0.4470 0.7410])
            plot3([DSD.x_v DSD.x_v], [-DSD.r_v DSD.r_v] , 3*[1 1], '*' , 'Color', [0.8500 0.3250 0.0980])
            
            % Plot Velocity Magnitude Field
            surf(DSD.x_mesh, DSD.y_mesh, DSD.u_norm_mesh) ; shading flat;
            
            % Add Streamlines (Complicated, because we generate them in 3d to place them where we like!)
            z_streamlines = 3;
            n_streamlines = 30;
            streamline_options = [0.1 20000];                            % Step (0.1=default) and Max Points (10000=default, but insufficient for us!)
            
            starty = linspace(min(min(DSD.y_mesh)), max(max(DSD.y_mesh)), n_streamlines);
            startx = min(min(DSD.x_mesh)) * ones(1, n_streamlines);
            startz = z_streamlines * ones(size(starty));
            
            x_mesh_bis = zeros(size(DSD.x_mesh, 1) , size(DSD.x_mesh, 2) , 2); x_mesh_bis(:,:, 1) = DSD.x_mesh; x_mesh_bis(:,:, 2) = DSD.x_mesh;
            y_mesh_bis = zeros(size(DSD.y_mesh, 1) , size(DSD.y_mesh, 2) , 2); y_mesh_bis(:,:, 1) = DSD.y_mesh; y_mesh_bis(:,:, 2) = DSD.y_mesh;
            u_mesh_bis = zeros(size(DSD.u_mesh, 1) , size(DSD.u_mesh, 2) , 2); u_mesh_bis(:,:, 1) = DSD.u_mesh; u_mesh_bis(:,:, 2) = DSD.u_mesh;
            v_mesh_bis = zeros(size(DSD.v_mesh, 1) , size(DSD.v_mesh, 2) , 2); v_mesh_bis(:,:, 1) = DSD.v_mesh; v_mesh_bis(:,:, 2) = DSD.v_mesh;
            
            z_mesh_bis =  ones(size(x_mesh_bis)) * z_streamlines; z_mesh_bis(:,:,1) = 0.9 * z_mesh_bis(:,:,1); z_mesh_bis(:,:,2) = 1.1 * z_mesh_bis(:,:,2);
            w_mesh_bis = zeros(size(x_mesh_bis));
            
            streamline(x_mesh_bis,y_mesh_bis, z_mesh_bis ,u_mesh_bis,v_mesh_bis, w_mesh_bis, startx,starty, startz, streamline_options);
            
            view(2); % axis([x_min x_max y_min y_max])
            colorbar; caxis([0 2]);
            
            % Decorate
            xlabel('x'); ylabel('y');
            title(['Adaptive Wake Solution - RES_{\Gamma} = ' , num2str(DSD.RES_stretch) , ...
                                         ' - RES_{shape} = '  , num2str(DSD.RES_shape  ) ]);
            legend(['Actuator (f_1^{av}=' , num2str(DSD.f1_average , 3) , ' , d_1=' , num2str(DSD.d1 , 3) , ')'] , ...
                   ['Vortices (x_v=' , num2str(DSD.x_v, 3) , ' , r_v=' , num2str(DSD.r_v, 3) , ')'] , ...
                    'Velocity Magnitude (|U|)'      , ...
                    'Velocity Field Streamlines'   , ...
                    'Location', 'SouthWest');
        end 
    end
    
end

