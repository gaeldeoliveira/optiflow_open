classdef dogoro_double_actuator_solver < handle
    %DOGORO_DOUBLE_ACTUATOR_SOLVER Manager Class in charge of initializing,
    %   running and postprocessing a double actuator (constant loading)
    %   inviscid incompressible flow case
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                     %
    %   Dogoro -    A simple specialized 2d Vorticity Equation Solver for %
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
        % Free-Stream Charactheristics
        u_inf               = 1;                            % Unperturbed Wind Speed
        rho                 = 1;                            % Density

        % Rotor Inputs
        f1                                                  % Upstream   Actuator Force Density
        f2                                                  % Downstream Actuator Force Density
        Dx                  = 0.525;                        % Distance between actuators (1 for good looking config, 40 to tend towards optimal solution)        
        d1                  = 1;                            % 1st Actuator Diameter
                        
        % Parametrization Inputs (for double actuator streamtube)
        theta_start         = 0;                            % Zero means wake starts at actuator
        theta_end           = 0.9999*pi()/2;                % pi()/2 means wake extends down to infinity
        n_stances           = 80;                           % n_stances = n_panels + 1 (usually 20 to 30 more than the number of nearwake stances!)
        n_farwake_stances   = 0;                            % Add smooth far wake deformations
        cos_exp             = 4;                            % Lumping factor for (super)cosine discretization (2 corresponds to cosine, 3-5 are more lumped near actuator!)
        n_nearwake_stances  = 30;                           % Between 30 and 4 here        
        alpha_start         = 0;                            % Zero means wake starts at actuator
        alpha_end           = pi();                         % pi()/2 means wake extends down to infinity
        beta_vector_guess   = [1 1];                        % Coefficient vector (ones(...) means straight wake)
        xs                  = 0;                            % Origin of h-theta system in x-y plane
        ys                  = 0;                            %

        % Numerical Parameters
        N_iterations        = 200000;                       % Maximum Number of Iterations
        information_period  = 10000;                        % Display information/plot every information_period iterations
        max_RES_stretch     = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
        max_RES_shape       = 1e-10;                        % Maximum Shape residual (volume flow accross segments)
        % Numerical Parameters 
        relax_stretch       = 0.001;                        % Relaxation Factor for the Extension Residual
        relax_shape         = 0.001;                        % Relaxation Factor for the Shape Residual
        % Computational Parameters
        use_parallel        = false;                         % If true, compute induction in parallel
        
        % % Processed Parameters
        r1                                                  % Radius of  Upstream Actuator
        S1                                                  % Surface of Upstream Actuator
        l_nearwake                                          % Lenght of nearwake
        n_panels                                            % Number of Panels on Mainwake
        
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
        
        % Global Residual Values
        RES_stretch                                         % Wake vorticity strenght residual 
        RES_shape                                           % Wake shape misalignment residual
        
        % % Postprocessing
        % Machine Parameters
        x1                                                  % x Position        (First  Actuator )
        a1                                                  % Induction Factor  (Absolute)
        m1                                                  % Mass Flux         (First  Actuator )
        p1                                                  % Power             (First  Actuator )
        Cp1                                                 % Power Coefficient (First  Actuator )
        % Downstream Actuator
        x2                                                  % x Position        (Second Actuator )
        r2                                                  % Radius            (Second Actuator )
        d2                                                  % Diameter          (Second Actuator )
        S2                                                  % Surface           (Second Actuator )
        a2                                                  % Induction Factor  (Absolute)
        m2                                                  % Mass Flux         (Second Actuator )
        p2                                                  % Power             (Second Actuator )
        Cp2                                                 % Power Coefficient (Second Actuator )
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
        
        % % Results Bundle
        RES                                                 % Bundles key machine parameters for compact storage to explore performance trends
        
    end
    
    methods
        function DSD = dogoro_double_actuator_solver(f1, f2, Dx)
            % Constructor Method for the Double Actuator Solver Class
            % % Import Actuator Parameters
            % First Actuator Force Density
            DSD.f1 = f1;
            % Second Actuator Force Density
            DSD.f2 = f2;
            % Distance between Actuators
            DSD.Dx = Dx;
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
            % And surface (2d surface is just diameter/width)
            DSD.S1 = DSD.d1;
            % Lenght of nearwake
            DSD.l_nearwake          = DSD.Dx;
        end
        
        function define_discretization_and_initial_geometry(DSD)
            % % Generate geometry and stance steps
            % Addapt
            beta_vector = DSD.r1*DSD.beta_vector_guess;                          % Scale Coefficients to radius
            
            % % Generate range of theta values for panel definition
            % (almost always between 0 and pi()/2)
            theta_range = linspace(DSD.theta_start, DSD.theta_end, ...
                                   DSD.n_stances-DSD.n_farwake_stances-DSD.n_nearwake_stances);  
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
            
            % % Add inter-actuator discretization for double actuator
            %   streamtube (cosine lumping!)
            if DSD.n_nearwake_stances > 0
                alpha_range = linspace(DSD.alpha_start, DSD.alpha_end,  DSD.n_nearwake_stances+1);  % Range of Theta Values (always between 0 and pi())
                %cos_discretization = (cos(alpha_range).^2).^(1/IM.cos_exp).*sign(cos(alpha_range));
                x_nearwake  = - 0.5 * (cos(alpha_range) + 1) * DSD.l_nearwake;
                x_nearwake  = x_nearwake(1:end-1); % (exclude last point at 0!)
                % y_nearwake = r2 * ones(size(x_nearwake));
                x = [x_nearwake(:) ; x(:)];
                y = DSD.r1 * ones(size(x));
            end
                                    
            % % Generate discretization from initial geometry/stance steps
            % Determine number of panels
            DSD.n_panels = DSD.n_stances - 1;
            
            % Now segment stances into panels (this could also be done earlier)
            DSD.x_start = x(1:(end-1));
            DSD.y_start = y(1:(end-1));
            
            DSD.x_end   = x(2:end);
            DSD.y_end   = y(2:end);
            
            % And duplicate into upper and lower wake filament
            DSD.x_start = [DSD.x_start(:) ;   DSD.x_start(:)];
            DSD.y_start = [DSD.y_start(:) ; - DSD.y_start(:)];
            
            DSD.x_end   = [DSD.x_end(:)   ;   DSD.x_end(:)];
            DSD.y_end   = [DSD.y_end(:)   ; - DSD.y_end(:)];
            
        end
        
        function define_reference_panel_strenghts(DSD)
            % % % Initial guess                        
            
            % % Generate array of forces
            f_panels = ones(1, DSD.n_panels);
            % Make early wake, due to 1st actuator only
            f_panels(1:DSD.n_nearwake_stances) = DSD.f1;
            % Make late wake, superposition of 1st and 2nd actuator
            f_panels(DSD.n_nearwake_stances+1:DSD.n_panels) = DSD.f1 + DSD.f2;
            
            % % Generate array of reference circulation densities
            % The reference circulation density of the panels is based on
            % the situation in which the tangential speed on the panel
            % equals free-stream.
            %
            % The equal circulation density is found by the solver, which
            % scaling this reference until the extension (stretch) residual
            % is minimized 
            % 
            % Upper side:
            DSD.gamma0(1:DSD.n_panels)       = f_panels./DSD.u_inf;
            % Lower side:                         (symmetric of upper side)
            DSD.gamma0((DSD.n_panels+1):(2*DSD.n_panels)) = ...
                                            - DSD.gamma0(1:DSD.n_panels);
            % Turn into Column Vector
            DSD.gamma0 = DSD.gamma0(:);
                        
        end
        
        function formulate_initial_guesses(DSD)
            % Set initial guess as current guess
            DSD.gamma = DSD.gamma0;
            % Allocate vector of normal velocities
            DSD.t_n = zeros(size(DSD.x_end));
        end
        
        function create_vortex_element_object(DSD)
            if not(DSD.use_parallel)
                DSD.VS    = constant_strenght_vortex_segment_2d(DSD.x_start, DSD.y_start, DSD.x_end, DSD.y_end, DSD.gamma);
            else
                DSD.VS    = constant_strenght_vortex_segment_2d_parallel(DSD.x_start, DSD.y_start, DSD.x_end, DSD.y_end, DSD.gamma);
            end
        end
        
        function run_solver(DSD)
            figure(9); hold on; grid on;
            for n=1:DSD.N_iterations
                % % Compute Movement on Panel Centers
                %   Zero-out x movement, this makes
                %   convergence a lot slower but ensures that we preserve a very fine
                %   discretization in critical areas.
                dx_center = 0 * DSD.t_n.* DSD.VS.x_n_unit_vector * DSD.relax_shape;
                dy_center =     DSD.t_n.* DSD.VS.y_n_unit_vector * DSD.relax_shape;
                
                %   Bypass deformation of last (most downstream) panel to prevent
                %   wake from closing on itself due to starting vortex, this is an early
                %   approach, can be refined with semi-infinite vortex segment later on!
                dx_center(DSD.n_panels) = 0;    dx_center(end) = 0;
                dy_center(DSD.n_panels) = 0;    dy_center(end) = 0;
                
                % % Move Panel Ends
                %   Transfer from center to ends (correct from first order
                %   point of view). Movement accumulates as we go forward. 
                %   Upper side first:
                DSD.x_end(1:DSD.n_panels) = DSD.x_end(1:DSD.n_panels) + cumsum(dx_center(1:DSD.n_panels));
                DSD.y_end(1:DSD.n_panels) = DSD.y_end(1:DSD.n_panels) + cumsum(dy_center(1:DSD.n_panels));
                %   Lower side second:
                DSD.x_end((DSD.n_panels+1):end) = DSD.x_end((DSD.n_panels+1):end) + cumsum(dx_center((DSD.n_panels+1):end));
                DSD.y_end((DSD.n_panels+1):end) = DSD.y_end((DSD.n_panels+1):end) + cumsum(dy_center((DSD.n_panels+1):end));
                
                % % Synchronize starts with ends
                %   Make x-starts meet x-ends
                DSD.x_start( 2:DSD.n_panels)        =  DSD.x_end(               1:(DSD.n_panels-1));
                DSD.x_start((DSD.n_panels+2):end)   =  DSD.x_end((DSD.n_panels+1):(end-1));
                %   Make y-starts meet y-ends
                DSD.y_start( 2:DSD.n_panels)        =  DSD.y_end(               1:(DSD.n_panels-1));
                DSD.y_start((DSD.n_panels+2):end)   =  DSD.y_end((DSD.n_panels+1):(end-1));
                
                % % Update VS object with new wake shape definition and panel strenghts
                DSD.VS.update_independent_fields(DSD.x_start, DSD.y_start, DSD.x_end, DSD.y_end, DSD.gamma);
                
                % % Compute Induced Speeds
                %   Speeds induced by wake, on panel centerpoints
                [~ , u_ctr , v_ctr ] = induced_speed_on_many_points(DSD.VS, DSD.VS.x_center, DSD.VS.y_center); %Evaluate speed at the center of the first segment (valid for multi-segment, with high resolution)
                %    Tangential Speed Component
                DSD.s_l = DSD.u_inf .* DSD.VS.x_l_unit_vector + (u_ctr .* DSD.VS.x_l_unit_vector + v_ctr .* DSD.VS.y_l_unit_vector);
                %     Normal Speed Component
                DSD.t_n = DSD.u_inf .* DSD.VS.x_n_unit_vector + (u_ctr .* DSD.VS.x_n_unit_vector + v_ctr .* DSD.VS.y_n_unit_vector);
                
                % Compute Residuals
                t_n_RES     = DSD.t_n( [1:(DSD.n_panels-1) , (DSD.n_panels+1):(2*DSD.n_panels-1)]);      % Exclude last panels!
                l_RES       = DSD.VS.l([1:(DSD.n_panels-1) , (DSD.n_panels+1):(2*DSD.n_panels-1)]);
                
                % % Compute Extension (stretch) Residual and Update Circulation Density
                %   Find Local Extension residual and its derivative
                RESl_stretch        =   (DSD.gamma - DSD.gamma0./DSD.s_l).^2 .* DSD.VS.l;
                RESl_stretch_dgamma = 2*(DSD.gamma - DSD.gamma0./DSD.s_l)    .* DSD.VS.l;
                %   Update Circulation Density
                DSD.gamma           = DSD.gamma - DSD.relax_stretch * RESl_stretch ./ RESl_stretch_dgamma;
                % Add clause to set NaNs to zero to deal with cases in
                % RESl_stretch_dgamma is exactly 0 (in case we have panels
                % without strenght! this does not enforce anything
                % definitive, the solver can change it again)
                DSD.gamma(isnan(DSD.gamma)) = 0;
                %     % Optional (for high loading)
                %     gamma(n_panels)     = gamma(n_panels-1);
                %     gamma(end)          = gamma(end-1);
                
                % % Compute Global Residuals
                DSD.RES_shape           = sum(t_n_RES.^2 .* l_RES);
                DSD.RES_stretch         = sum(RESl_stretch);
                
                % Display iteration information if demanded
                if or(n==1, round(n/DSD.information_period) == n/DSD.information_period)
                    DSD.display_periodic_information(n);
                end
                
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
            plot(n, log10(DSD.RES_shape), 'xk', n, log10(DSD.RES_stretch), 'ok');
            xlabel('N_{iter}'); ylabel('Residuals'); grid on;
            legend('Res_{shape}', 'Res_{\gamma}' );
            drawnow;
        end
        
        function compute_machine_parameters(DSD)            
            % Find Position and Radius of Second Actuator
            DSD.x2 = DSD.VS.x_start(DSD.n_nearwake_stances+1);
            DSD.r2 = DSD.VS.y_start(DSD.n_nearwake_stances+1);
            DSD.d2 = DSD.r2*2;
            DSD.S2 = DSD.d2;
            % Compute induction, mass flux and power on surface of second actuator
            [DSD.a2 , DSD.m2, DSD.p2] = DSD.compute_induction_on_cross_section(DSD.x2, DSD.r2, DSD.f2);
            % And finally reach power coefficient!
            DSD.Cp2  = DSD.p2 / (0.5*DSD.rho*DSD.S2*(DSD.u_inf^3));
            
            % Find Position and Radius of First Actuator
            DSD.x1 = DSD.VS.x_start(1);
            DSD.r1 = DSD.VS.y_start(1);
            DSD.d1 = DSD.r1*2;
            DSD.S1 = DSD.d1;
            % Compute induction, mass flux and power on surface of second actuator
            [DSD.a1 , DSD.m1, DSD.p1] = DSD.compute_induction_on_cross_section(DSD.x1, DSD.r1, DSD.f1);
            % And finally reach power coefficient!
            DSD.Cp1  = DSD.p1 / (0.5*DSD.rho*DSD.S2*(DSD.u_inf^3));
            
            % Now Compute Total Power and Complete Power Coefficient
            DSD.p    = DSD.p1 + DSD.p2;            
            DSD.Cp   = DSD.p  / (0.5*DSD.rho*DSD.S2*(DSD.u_inf^3));
            
            
            %% Display Diagnostic Text
            % Inputs
            disp('-- Problem Inputs --')
            disp(['d_1      = ' , num2str(DSD.d1) ]);
            disp(['f_1      = ' , num2str(DSD.f1) ]);
            disp(['f_2      = ' , num2str(DSD.f2) ]);
            disp(['D_x      = ' , num2str(DSD.Dx) ]);
            % Force
            disp('-- Actuator Size and Power --')
            disp(['r_1      = ' , num2str(DSD.r1) ]);
            disp(['r_2      = ' , num2str(DSD.r2) ]);
            disp(['p_1      = ' , num2str(DSD.p1) ]);
            disp(['p_2      = ' , num2str(DSD.p2) ]);
            
            % Power Coefficients
            disp('-- Power Coefficients --')
            disp(['Cp_1     = ' , num2str(DSD.Cp1) ]);
            disp(['Cp_2     = ' , num2str(DSD.Cp2) ]);
            disp(['Cp_total = ' , num2str(DSD.Cp ) ]);
            % Speed on Actuator
            disp('-- Average Normal Speed on Actuators --')
            disp(['u_1      = ' , num2str(DSD.m1/DSD.rho) ]);
            disp(['u_2      = ' , num2str(DSD.m2/DSD.rho)]);
            % Solution Parameters
            disp('-- Solution Residuals --')
            disp(['Strecht  = ' , num2str(DSD.RES_stretch) ]);
            disp(['Shape    = ' , num2str(DSD.RES_shape)   ]);
            
            %% Save to Results Bundle
            DSD.RES.d1         = DSD.d1;
            DSD.RES.f1         = DSD.f1;
            DSD.RES.f2         = DSD.f2;
            DSD.RES.Dx         = DSD.Dx;
            
            DSD.RES.r1         = DSD.r1;
            DSD.RES.r2         = DSD.r2;
            DSD.RES.p1         = DSD.p1;
            DSD.RES.p2         = DSD.p2;
            
            DSD.RES.Cp1        = DSD.Cp1;
            DSD.RES.Cp2        = DSD.Cp2;
            DSD.RES.Cp         = DSD.Cp ;
            
            DSD.RES.u1         = DSD.m1/DSD.rho;
            DSD.RES.u2         = DSD.m1/DSD.rho;
            DSD.RES.a1         = DSD.a1;
            DSD.RES.a2         = DSD.a2;
            
            DSD.RES.RES_stretch = DSD.RES_stretch;
            DSD.RES.RES_shape   = DSD.RES_shape;
            DSD.RES.VS          = DSD.VS;
            
            
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
            plot3([DSD.x2 DSD.x2], [-DSD.r2 DSD.r2] , 3*[1 1], '-' , 'Color', [0.8500 0.3250 0.0980])
            
            % Plot Velocity Magnitude Field
            surf(DSD.x_mesh, DSD.y_mesh, DSD.u_norm_mesh) ; shading interp;
            
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
            
            streamline(x_mesh_bis,y_mesh_bis, z_mesh_bis ,u_mesh_bis,v_mesh_bis, w_mesh_bis, startx,starty, startz, streamline_options)            
            
            view(2); % axis([x_min x_max y_min y_max])
            colorbar; % caxis([0 2]);
            
            % Decorate
            xlabel('x'); ylabel('y');
            title(['Adaptive Wake Solution - RES_{\Gamma} = ' , num2str(DSD.RES_stretch) , ...
                                         ' - RES_{shape} = '  , num2str(DSD.RES_shape  ) ]);
            legend(['Actuator 1  (f_1=' , num2str(DSD.f1, 3) , ' , d_1=' , num2str(DSD.d1, 3) , ')'] , ...
                   ['Actuator 2  (f_2=' , num2str(DSD.f2, 3) , ' , d_2=' , num2str(DSD.d2, 3) , ')'] , ...
                    'Velocity Magnitude (|U|)'      , ...
                    'Velocity Field Streamlines'   , ...
                    'Location', 'SouthWest');
        end
        
    end
    
end

