classdef experimental_data_loader2 < handle
    %EXPERIMENTAL DATA LOADER helps us load and use Ricardo's processed 
    % experimental data 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Integral Boundary Layer Integrator (ODE)
    %           Plasma Development Tool
    %
    %       August 2014, GNU-GPLv3 or later
    %       Gael de Oliveira, Ricardo Pereira
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %       Load Pereira's Processed Experimental Data
    %               following:
    %                   some reference here
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % Case Properties
        U_inf       = 20                                % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)
        rho_inf     = 1.225                             % [kg/m3]   Free Stream Density
        nu_inf      = 1.461e-5                          % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
        L           = 1                                 % [m]       Lenght Scale
        
        % Actuator Specific Definitions
        T_p         = .0025                             % [m]       Single Actuator DBD Force Field Thickness
        L_p         = .0175                             % [m]       Single Actuator DBD Force Field Lenght
        CFT         = -.15/3                             % [N/m]     Single Actuator DBD Total Force per unit span
        
        N_p         = 3                                 % [.]       Number of Actuators
        D_p         = .0400                             % [m]       Spacing between actuator Leading Edges


        % Raw Data [SI units]
        X                                               % [m]       Body attached tangential coordinate (usually known as S)
        Ue                                              % [m/s]     Edge Velocity  
        delta_1                                         % [m]       Displacement thickness
        delta_2                                         % [m]       Momenum thickness
        V_BL                                            % [m]       Velocity vectors in BL
        X0_p0                                           % [m]       (?) Start of first plasma actuator in body attached tangential coordinate 
        
        % Processed data [Adimensional]
        xi_vector                                       % [--]      Adimensional tangential coordinate (xi = X / L)
        ue_vector                                       % [--]      Adimensional edge speed (ue = Ue / U_inf)
        dstr_vector                                     % [--]      Adimensional Displacement Thickness (dstr = delta_1 / L)
        theta_vector                                    % [--]      Adimensional Momentum Thickness (theta = delta_2 / L)
        u_inf_over_nu_inf                               % [m^(-1)]  Viscosity Scale (Numerically equal to the unit Lenght Scale Reynolds)
        
        % Generated Plasma Descriptor Interfaces        
        X0_p_vector                                     % [m]       Vector of DBD Actuator Force Field Starting Stances (one entry per actuator!)
        L_p_vector                                      % [m]       Vector of DBD Actuator Force Field Lenghts (one entry per actuator!)
        T_p_vector                                      % [m]       Vector of DBD Actuator Force Field Thicknesses (one entry per actuator!)
        CFT_p_vector                                    % [m]       Vector of DBD Actuator Force Field Thicknesses (one entry per actuator!)
        
        % IO Support
        experiment_folder_path = '../Processing_PIV/'   % Path to experiment folder
        experiment_file                                 % Experiment file Name
        experiment_data                                 % Experiment data structure        
    end
    
    methods
        function EDL = experimental_data_loader2()
            % Constructor Method (explicit declaration! just because we are onanists)
        end
        
        function load_experiment_data(EDL, experiment_file)
            % function load_experiment_data(EDL, experiment_file)
            %   Handles file IO operations to load experiment data from MAT
            %   file
            
            % Load experimental data variables of interest into a structure
            EDL.experiment_file = experiment_file;
            EDL.experiment_data = load([EDL.experiment_folder_path experiment_file]);
            % 'airfoil_ref' , 'U0' , 'dstar' , 'theta', 'x_DBD_airfoil' , 'y_DBD_airfoil ' , 'V_BL');
            
            % Extract Experiment Data
            EDL.extract_experiment_data();
            % Process Experiment Data
            EDL.process_experiment_data();
            % Generate Plasma Descriptor Interfaces
            EDL.generate_plasma_descriptor_interfaces();
            
        end
        
        function extract_experiment_data(EDL)
            % function extract_experiment_data(EDL)
            %   Extracts experiment results from the experiment_data
            %   structure loaded through the load_experiment_data function
            
            % Store edge velocity as is
            EDL.Ue = EDL.experiment_data.U0;
            
            % Store thicknesses as they come
            EDL.delta_1 = EDL.experiment_data.dstar;
            EDL.delta_2 = EDL.experiment_data.theta;
            
            % Now process the airfoil coordinates 
            %       From (x,y) to -> (s, n)
            %            [mm]  to [m]
            % (CAD files store in milimeters)
            % Get airfoil coordinates of interest
            x = EDL.experiment_data.airfoil_ref(:,1) / 1000;                       % Airfoil x coordinates, from [mm] to [m]
            y = EDL.experiment_data.airfoil_ref(:,2) / 1000;                       % Airfoil y coordinates, from [mm] to [m]
            % Compute Steps between coordinates
            dx = x(2:end) - x(1:(end-1));
            dy = y(2:end) - y(1:(end-1));            
            % Compute Norm of steps
            ds = sqrt(dx.^2 + dy.^2);
            % Set right array format (robust)
            ds = ds(:)';
            % Generate s coordinate data (tangential to body)
            s = [0 cumsum(ds)];            
            % Store
            EDL.X = s;
            
            % % Find the Position of the beginning of the first plasma actuator
            % Look for minimum distance to coordinates!
            [~ , i_d_min] = min((EDL.experiment_data.x_DBD_airfoil - EDL.experiment_data.airfoil_ref(:,1)).^2 + ...
                                (EDL.experiment_data.y_DBD_airfoil - EDL.experiment_data.airfoil_ref(:,2)).^2);
            % Store position!
            EDL.X0_p0 = s(i_d_min);
            
            % % Store Velocity Profiles for reference
            EDL.V_BL = EDL.experiment_data.V_BL;
            
        end
        
        function process_experiment_data(EDL)
            % Generate the Processed data fields, from the Raw Data Fields
            % and the Case Definition Variables
            
            % [--]      Adimensional tangential coordinate (xi = X / L)
            EDL.xi_vector           = EDL.X ./ EDL.L;
            
            % [--]      Adimensional edge speed (ue = Ue / U_inf)
            EDL.ue_vector           = EDL.Ue ./ EDL.U_inf;
            
            % [--]      Adimensional Displacement Thickness (dstr = delta_1 / L)
            EDL.dstr_vector         = EDL.delta_1 ./ EDL.L;
            
            % [--]      Adimensional Momentum Thickness (theta = delta_2 / L)
            EDL.theta_vector        = EDL.delta_2 ./ EDL.L;
            
            % Viscosity Scale
            EDL.u_inf_over_nu_inf   = EDL.U_inf ./ EDL.nu_inf;                        
        end
        
        function generate_plasma_descriptor_interfaces(EDL)
            
            % Allocate Plasma Descriptor Interface Arrays
            EDL.X0_p_vector     = zeros(1,EDL.N_p);
            EDL.L_p_vector      = zeros(1,EDL.N_p);
            EDL.T_p_vector      = zeros(1,EDL.N_p);
            EDL.CFT_p_vector    = zeros(1,EDL.N_p);
            
            % Fill Plasma Descriptor Interface Arrays in
            for n_p = 1:EDL.N_p
                % Start of Each plasma actuator, depends on spacing and
                % plasma actuator number
                EDL.X0_p_vector(n_p)    = EDL.X0_p0 + (n_p -1) * EDL.D_p;
                % Generate Repetitive plasma descriptor!
                EDL.L_p_vector(n_p)     = EDL.L_p;
                EDL.T_p_vector(n_p)     = EDL.T_p;
                EDL.CFT_p_vector(n_p)   = EDL.CFT;                            
            end
                                               
        end
        
        function plot(EDL)
            fig = figure();
            hax=axes; 
            plot(EDL.X, EDL.delta_1, EDL.X, EDL.delta_2)
            grid on
            xlabel('X [m] - Body Attached Tangential Coordinate')
            ylabel('Thicknesses [m]')
            legend('delta_1 [m] - Displacement Thickness' , 'delta_2 [m] - Momentum Thickness')
            
            line([EDL.X0_p EDL.X0_p],get(hax,'YLim'),'Color',[1 0 0])
        end

    end
    
end

