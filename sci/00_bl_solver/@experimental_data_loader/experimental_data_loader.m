classdef experimental_data_loader < handle
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
        experiment_folder_path = '../Processing_PIV/'   % Path to experiment folder
        experiment_file                                 % Experiment file Name
        experiment_data                                 % Experiment data structure  

        X                                               % [m]   Body attached tangential coordinate (usually known as S)
        Ue                                              % [m/s] Edge Velocity  
        delta_1                                         % [m]   displacement thickness
        delta_2                                         % [m]   momenum thickness
        V_BL                                            % [m]   velocitie vectors in BL
        X0_p                                            % [m]   Start(?) of plasma in body attached tangential coordinate 
        
        

        
    end
    
    methods
        function EDL = experimental_data_loader()
            
        end
        
        function load_experiment_data(EDL, experiment_file)
            % Load experimental data variables of interest into a structure
            EDL.experiment_file = experiment_file;
            EDL.experiment_data = load([EDL.experiment_folder_path experiment_file]);
            % 'airfoil_ref' , 'U0' , 'dstar' , 'theta', 'x_DBD_airfoil' , 'y_DBD_airfoil ' , 'V_BL');
            
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
            
            % % Find the Position of the plasma actuator
            % Look for minimum distance to coordinates!
            [~ , i_d_min] = min((EDL.experiment_data.x_DBD_airfoil - EDL.experiment_data.airfoil_ref(:,1)).^2 + ...
                                (EDL.experiment_data.y_DBD_airfoil - EDL.experiment_data.airfoil_ref(:,2)).^2);
            % Store position!
            EDL.X0_p = s(i_d_min);
            
            % % Store Velocity Profiles for reference
            EDL.V_BL = EDL.experiment_data.V_BL;
            
        end
        
        function process_experiment_data(EDL, experiment_file)
            
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

