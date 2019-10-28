classdef experimental_case < handle
    %EXPERIMENTAL CASE is a handle class for storing the results and inputs
    % of a wind tunnel experiment in which the lift, drag and pitch moment
    % coefficients of an airfoil have been measured at constant speed, mach
    % and Reynolds number. 
    %
    % The usage of handle is favored to limit data exchange needs, but this
    % may be waived at a later point if there are issues with parallelism.
    % 
    %class 
    %   It is expected that lift data always exists, while drag and
    %   pitch moment coefficient data may or may not be present. Each
    %   output type 
    %
    %   In the future, this class will be extended to cover other types of
    %   experiments
    %
    %   Precautions for unusual configurations will be inserted at a later
    %   stage
    
    properties
        % For templating
        uuid                        % [uuid  ] - UUID of experimental case 
        EDB                         % [class ] - Experimental case database to which this case belongs
        EDB_uuid                    % [uuid  ] - UUID of experimental case database
        
        % Specific: Replicated data (for UI search and display purposes only)
        airfoil_name                % [string] - String with friendly airfoil name (for UI search and display purposes only)
        data_source                 % [string] - String with source of data (for UI search and display purposes only)
        Re                          % [adim.] - Reynolds Number (for UI search and display purposes only)
    end
    
    properties (SetAccess= private)
        % For templating
        IDkind = 'airfoil_polar'    % [string] - Protected string describes the type of case this describes
        
        % Specific: Structure data stores
        airfoil_description         % [obj.  ] - Description (templated) of airfoil shape and data source
        polar_conditions            % [obj.  ] - Description (templated) of aerodynamic polar conditions
        processed_polars            % [obj.  ] - Description (templated) of aerodynamic polar data
        
        % Processed Fields and Objects (may enter future template)
        flag_alcl_polar_present     % [bool. ] - Indicates (when true) that the experimental case includes data for cl vs alpha
        flag_clcd_polar_present     % [bool. ] - Indicates (when true) that the experimental case includes data for cl vs cd
        flag_alclcd_polar_present   % [bool. ] - Indicates (when true) that the experimental case includes data for cl and cd vs alpha
        flag_alcm_polar_present     % [bool. ] - Indicates (when true) that the experimental case includes data for cm vs alpha
        
        
        % Standardized polar data (corrected for alpha_offset and chord_ratio)
        al_alcl_polar_std           % [deg   ] - Column array with standardized angles of attack at which lift was measured
        cl_alcl_polar_std           % [adim. ] - Column array with standardized lift coefficients
        
        cl_clcd_polar_std           % [adim. ] - Column array with standardized lift coefficients
        cd_clcd_polar_std           % [adim. ] - Column array with standardized drag coefficients
        
        al_alclcd_polar_std         % [deg.  ] - Column array with standardized angles of attack
        cl_alclcd_polar_std         % [adim. ] - Column array with standardized lift coefficients
        cd_alclcd_polar_std         % [adim. ] - Column array with standardized drag coefficients
        
        al_alcm_polar_std           % [deg   ] - Column array with standardized angles of attack at which picht moment was measured
        cm_alcm_polar_std           % [adim. ] - Column array with standardized pitch moment coefficients
    end
    
    methods
        function EC = experimental_case(EDB)
            %Construct an instance of the experimental case class
            %   Instances of the experimental case class are constructed
            %   without providing data. Instead, they receive a hook to an
            %   an experimental case database object, which guides them on
            %   how to process data in a standardized way. 
            %
            %   Even though the EDB is not a parent in the formal sense of
            %   the world, it can be thought of as such because it governs
            %   the behavior of the EC to a large extent.
            
            % Store hook to EDB object 
            EC.EDB = EDB;
            EC.EDB_uuid = EC.EDB.uuid;
            % Make uuid of experimental case
            EC.uuid = EC.EDB.make_a_uuid();
            % Initialize appropriate flags to boolean 
            %EC.flag_complete_insert     = false;
            %EC.flag_complete_postproc   = false;
            EC.flag_alcl_polar_present   = false;
            EC.flag_clcd_polar_present   = false;
            EC.flag_alclcd_polar_present = false;
            EC.flag_alcm_polar_present   = false;
        end
        
        function load_airfoil_description(EC, airfoil_description)
            % Store Airfoil Description
            EC.airfoil_description = airfoil_description;
            % Repeat airfoil name (for display purposes only)
            EC.airfoil_name = airfoil_description.airfoil_name;
        end

        function load_polar_conditions(EC, polar_conditions)
            % Store Polar Conditions Description
            EC.polar_conditions = polar_conditions;
            % Repeat Reynolds number (for display purposes only)
            EC.Re = polar_conditions.Re;
        end

        function load_processed_polars(EC, processed_polars)
            % Store processed polars
            EC.processed_polars = processed_polars;
            % Set flags
            EC.flag_alcl_polar_present   = false;
            EC.flag_clcd_polar_present   = false;
            EC.flag_alclcd_polar_present = false;
            EC.flag_alcm_polar_present   = false;
            % Set presence flag for alpha, cl polar
            if and(not(isempty(processed_polars.al_alcl_polar)), not(isempty(processed_polars.cl_alcl_polar)))
                EC.flag_alcl_polar_present   = true;
            end
            % Set presence flag for alpha, cm polar
            if and(not(isempty(processed_polars.al_alcm_polar)), not(isempty(processed_polars.cm_alcm_polar)))
                EC.flag_alcm_polar_present   = true;
            end
            
            % Set presence flag for cl, cd polar
            if and(not(isempty(processed_polars.cl_clcd_polar)), not(isempty(processed_polars.cd_clcd_polar)))
                EC.flag_clcd_polar_present   = true;
            end
            % Set presence flag for alpha, cl, cd polar
            if 3 == (not(isempty(processed_polars.al_alclcd_polar)) + not(isempty(processed_polars.cl_alclcd_polar)) + not(isempty(processed_polars.cd_alclcd_polar)))
                EC.flag_alclcd_polar_present = true;
            end
            % Finally, get some metadata for display purposes only
            EC.data_source = EC.processed_polars.source;
        end

        function standardize_polars(EC)
            % This function standardizes polar data for the effect of
            % scaled and rotated airfoil coordinates 
            %
            % TODO: check that sign of alpha_offset is correct!
            %
            
            % First, store alpha_offset and chord_ratio
            alpha_offset = EC.airfoil_description.theta;
            scale_factor  = EC.airfoil_description.scale_factor;
            
            % Standardize Alpha, cl polar
            EC.al_alcl_polar_std   = EC.processed_polars.al_alcl_polar   + alpha_offset;
            EC.cl_alcl_polar_std   = EC.processed_polars.cl_alcl_polar   * scale_factor;
            % Standardize Alpha, cm polar
            EC.al_alcm_polar_std   = EC.processed_polars.al_alcm_polar   + alpha_offset;
            EC.cm_alcm_polar_std   = EC.processed_polars.cm_alcm_polar   * scale_factor;
            % Standardize Cl, cd polar
            EC.cl_clcd_polar_std   = EC.processed_polars.cl_clcd_polar   * scale_factor;
            EC.cd_clcd_polar_std   = EC.processed_polars.cd_clcd_polar   * scale_factor;
            % Standardize Alpha, cl, cd polar
            EC.al_alclcd_polar_std = EC.processed_polars.al_alclcd_polar + alpha_offset;
            EC.cl_alclcd_polar_std = EC.processed_polars.cl_alclcd_polar * scale_factor;
            EC.cd_alclcd_polar_std = EC.processed_polars.cd_alclcd_polar * scale_factor;
            
        end

        function cl = cl_alpha(EC, alpha)
            cl = interp1(EC.al_alcl_polar_std  , EC.cl_alcl_polar_std  , alpha);
        end

        function cm = cm_alpha(EC, alpha)
            cm = interp1(EC.al_alcm_polar_std  , EC.cm_alcm_polar_std  , alpha);
        end

        function cd = cd_given_cl(EC, cl)
            cd = interp1(EC.cl_clcd_polar_std  , EC.cd_clcd_polar_std  , cl   );
        end

        function cd = cd_alpha(EC, alpha)
            cd = interp1(EC.al_alclcd_polar_std, EC.cd_alclcd_polar_std, alpha);
        end
        
        function n_datapoints = count_datapoints(EC)
            % Count number of datapoints in experimental case
            n_datapoints = length(EC.cl_alcl_polar_std) + length(EC.cd_alclcd_polar_std) + length(EC.cm_alcm_polar_std);
        end
        
        function plot(EC)
            figure(1)
            EC.processed_polars 
            subplot(122)
            plot(EC.processed_polars.al_alcl_polar   , EC.processed_polars.cl_alcl_polar   , 'x-'); hold on;
            plot(EC.processed_polars.al_alclcd_polar , EC.processed_polars.cl_alclcd_polar , 'o-'); grid on;
            xlabel('\alpha (deg)'); ylabel('C_l');
            subplot(121)
            plot(EC.processed_polars.cd_clcd_polar   * 1e4, EC.processed_polars.cl_clcd_polar   , 'x-'); hold on;
            plot(EC.processed_polars.cd_alclcd_polar * 1e4, EC.processed_polars.cl_alclcd_polar , 'o-'); grid on;
            xlabel('C_d x 10^4') ; ylabel('C_l');
            legend('Sorted Data' , 'Reinterpolated Data', 'Location', 'East');
            
            figure(2)
            plot(EC.airfoil_description.tx_coordinates    , EC.airfoil_description.tz_coordinates          ); hold on;
            plot(EC.airfoil_description.tx_coordinates_raw, EC.airfoil_description.tz_coordinates_raw, '--.'); grid on;
            [tx, tz] = EC.EDB.SD.generate_coordinates(160 , EC.airfoil_description.x);
            plot(tx, tz, '-.'); grid on;
            xlabel('x/c'); ylabel('y/c'); axis equal;
            legend('Standardized Airfoil Shape' , 'Raw Airfoil Shape', 'Regenerated Airfoil Shape');
        end
        
        
    end
    
end

