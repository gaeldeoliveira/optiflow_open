classdef TR824_reader
    %TR824_READER Collects static methods for reading and structuring data
    %from digitized TR824 report.
    
    properties
    end
    
    methods
        function T824R = TR824_reader()
            % Dummy Constructor
        end
    end
    
    methods(Static)
        function [line_str, line_data] = read_dataset_data_line(fid)
            try
                line_str      = fgetl(fid);
                line_data = sscanf(line_str, '%f', 3);
            catch
                % Signal EOL with empty return!
                line_str  = [];
                line_data = [];
            end
        end
        
        function [dataset_number, dataset_data] = read_dataset_data(fid)
            % Dataset number
            line_dataset_number = fgetl(fid);
            dataset_number      = sscanf(line_dataset_number, '%f');
            
            dataset_data   = 0;
            n_dataset_line = 1;
            while true
                [~, line_data] = TR824_reader.read_dataset_data_line(fid);
                if length(line_data) == 2
                    dataset_data(n_dataset_line, 1) = line_data(1);  %#ok<AGROW>
                    dataset_data(n_dataset_line, 2) = line_data(2);  %#ok<AGROW>
                elseif (length(line_data) == 3)
                    dataset_data(n_dataset_line, 1) = line_data(1);  %#ok<AGROW>
                    dataset_data(n_dataset_line, 2) = line_data(2);  %#ok<AGROW>
                    dataset_data(n_dataset_line, 3) = line_data(3);  %#ok<AGROW>
                else
                    % Finished Dataset!
                    % Deal with empty datasets if necessary
                    if dataset_data == 0
                        % Return empty data
                        dataset_data = [];
                        % Skip extra line of empty datasets
                        fgetl(fid);
                    end
                    % Break whatever the reason (succes or emptyness) for being here is!
                    break;
                end
                n_dataset_line = n_dataset_line + 1;
            end
        end
        
        function [datasets] = read_datasets_from_file(filename)
           % Reads file from digitized NACA-TR824 and returns numeric data
           % grouped in a TR824_datasets structure (templated).
            
           % Emmit message stating that we are reading file
           lg.msg(TR824_reader, ['Reading Data from file : ' filename]);
           % Open the polar file for reading ()
           fid = fopen(filename, 'r');
           
           % Read first header line
           line_first_header               = fgetl(fid); %#ok<NASGU>
           % Read first dataset separator
           line_first_separator            = fgetl(fid); %#ok<NASGU>
           
           % % Start reading datasets
           % Read alpha vs cl (1-6) and cm (7-12) datasets
           [dataset_A01_number, dataset_A01_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A02_number, dataset_A02_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A03_number, dataset_A03_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A04_number, dataset_A04_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A05_number, dataset_A05_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A06_number, dataset_A06_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A07_number, dataset_A07_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A08_number, dataset_A08_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A09_number, dataset_A09_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A10_number, dataset_A10_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A11_number, dataset_A11_data] = TR824_reader.read_dataset_data(fid);
           [dataset_A12_number, dataset_A12_data] = TR824_reader.read_dataset_data(fid);
           
           % Read first header line
           line_second_header               = fgetl(fid); %#ok<NASGU>
           % Read first dataset separator
           line_second_separator            = fgetl(fid); %#ok<NASGU>
           % Read Cl, Cd (1-6) and Cl, Cm (7-12) datasets
           [dataset_B01_number, dataset_B01_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B02_number, dataset_B02_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B03_number, dataset_B03_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B04_number, dataset_B04_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B05_number, dataset_B05_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B06_number, dataset_B06_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B07_number, dataset_B07_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B08_number, dataset_B08_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B09_number, dataset_B09_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B10_number, dataset_B10_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B11_number, dataset_B11_data] = TR824_reader.read_dataset_data(fid);
           [dataset_B12_number, dataset_B12_data] = TR824_reader.read_dataset_data(fid);
           
           % Close file
           fclose(fid);
           
           % Make a consistency check on data set numbers
           if (dataset_A01_number ==  1); lg.msg(TR824_reader, 'PASS: dataset_A01 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A01 number'); end
           if (dataset_A02_number ==  2); lg.msg(TR824_reader, 'PASS: dataset_A02 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A02 number'); end
           if (dataset_A03_number ==  3); lg.msg(TR824_reader, 'PASS: dataset_A03 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A03 number'); end
           if (dataset_A04_number ==  4); lg.msg(TR824_reader, 'PASS: dataset_A04 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A04 number'); end
           if (dataset_A05_number ==  5); lg.msg(TR824_reader, 'PASS: dataset_A05 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A05 number'); end
           if (dataset_A06_number ==  6); lg.msg(TR824_reader, 'PASS: dataset_A06 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A06 number'); end
           if (dataset_A07_number ==  7); lg.msg(TR824_reader, 'PASS: dataset_A07 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A07 number'); end
           if (dataset_A08_number ==  8); lg.msg(TR824_reader, 'PASS: dataset_A08 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A08 number'); end
           if (dataset_A09_number ==  9); lg.msg(TR824_reader, 'PASS: dataset_A09 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A09 number'); end
           if (dataset_A10_number == 10); lg.msg(TR824_reader, 'PASS: dataset_A10 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A10 number'); end
           if (dataset_A11_number == 11); lg.msg(TR824_reader, 'PASS: dataset_A11 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A11 number'); end
           if (dataset_A12_number == 12); lg.msg(TR824_reader, 'PASS: dataset_A12 number'); else; lg.err(TR824_reader, 'FAIL: dataset_A12 number'); end
           
           
           if (dataset_B01_number ==  1); lg.msg(TR824_reader, 'PASS: dataset_B01 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B01 number'); end
           if (dataset_B02_number ==  2); lg.msg(TR824_reader, 'PASS: dataset_B02 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B02 number'); end
           if (dataset_B03_number ==  3); lg.msg(TR824_reader, 'PASS: dataset_B03 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B03 number'); end
           if (dataset_B04_number ==  4); lg.msg(TR824_reader, 'PASS: dataset_B04 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B04 number'); end
           if (dataset_B05_number ==  5); lg.msg(TR824_reader, 'PASS: dataset_B05 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B05 number'); end
           if (dataset_B06_number ==  6); lg.msg(TR824_reader, 'PASS: dataset_B06 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B06 number'); end
           if (dataset_B07_number ==  7); lg.msg(TR824_reader, 'PASS: dataset_B07 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B07 number'); end
           if (dataset_B08_number ==  8); lg.msg(TR824_reader, 'PASS: dataset_B08 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B08 number'); end
           if (dataset_B09_number ==  9); lg.msg(TR824_reader, 'PASS: dataset_B09 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B09 number'); end
           if (dataset_B10_number == 10); lg.msg(TR824_reader, 'PASS: dataset_B10 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B10 number'); end
           if (dataset_B11_number == 11); lg.msg(TR824_reader, 'PASS: dataset_B11 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B11 number'); end
           if (dataset_B12_number == 12); lg.msg(TR824_reader, 'PASS: dataset_B12 number'); else; lg.err(TR824_reader, 'FAIL: dataset_B12 number'); end
           
           % % Make object to store and return datasets
           datasets = TR824_datasets();
           % First A datasets
           datasets.A01_data = dataset_A01_data;    % alpha, Cl at Re=3e6, clean configuration
           datasets.A02_data = dataset_A02_data;    % alpha, Cl at Re=6e6, clean configuration
           datasets.A03_data = dataset_A03_data;    % alpha, Cl at Re=9e6, clean configuration
           datasets.A04_data = dataset_A04_data;    % alpha, Cl at Re=6e6, rough configuration
           datasets.A05_data = dataset_A05_data;    % alpha, Cl at Re=6e6, clean configuration, split flap
           datasets.A06_data = dataset_A06_data;    % alpha, Cl at Re=6e6, rough configuration, split flap
           datasets.A07_data = dataset_A07_data;    % alpha, Cm at Re=3e6, clean configuration, around c/4 point
           datasets.A08_data = dataset_A08_data;    % alpha, Cm at Re=6e6, clean configuration, around c/4 point
           datasets.A09_data = dataset_A09_data;    % alpha, Cm at Re=9e6, clean configuration, around c/4 point
           datasets.A10_data = dataset_A10_data;    % alpha, Cm at Re=6e6, rough configuration, around c/4 point
           datasets.A11_data = dataset_A11_data;    % alpha, Cm at Re=6e6, clean configuration, around c/4 point, split flap
           datasets.A12_data = dataset_A12_data;    % alpha, Cm at Re=6e6, rough configuration, around c/4 point, split flap
           % Then B datasets
           datasets.B01_data = dataset_B01_data;    % Cl, Cd at Re=3e6, clean configuration
           datasets.B02_data = dataset_B02_data;    % Cl, Cd at Re=6e6, clean configuration
           datasets.B03_data = dataset_B03_data;    % Cl, Cd at Re=9e6, clean configuration
           datasets.B04_data = dataset_B04_data;    % Cl, Cd at Re=6e6, rough configuration
           datasets.B05_data = dataset_B05_data;    % Cl, Cd at Re=6e6, clean configuration, split flap
           datasets.B06_data = dataset_B06_data;    % Cl, Cd at Re=6e6, rough configuration, split flap
           datasets.B07_data = dataset_B07_data;    % Cl, Cm at Re=3e6, clean configuration, around mac point
           datasets.B08_data = dataset_B08_data;    % Cl, Cm at Re=6e6, clean configuration, around mac point
           datasets.B09_data = dataset_B09_data;    % Cl, Cm at Re=9e6, clean configuration, around mac point
           datasets.B10_data = dataset_B10_data;    % Cl, Cm at Re=6e6, rough configuration, around mac point
           datasets.B11_data = dataset_B11_data;    % Cl, Cm at Re=6e6, clean configuration, around mac point, split flap
           datasets.B12_data = dataset_B12_data;    % Cl, Cm at Re=6e6, rough configuration, around mac point, split flap
           % Metatada
           datasets.filename = filename;            % Name of filename from which TR824 data was read
        end
        
        function [processed_polars] = process_polar_triplet(alcl_polar_dataset, alcm_polar_dataset, clcd_polar_dataset)
            % Returns a structure, processed_polars, with fields:
            % First : alpha, cl polar
            %processed_polars.al_alcl_polar    % Angle of attack (deg)
            %processed_polars.cl_alcl_polar    % Lift coefficient
            % Second: cl, cd polar
            %processed_polars.cl_clcd_polar    % Lift coefficient
            %processed_polars.cd_clcd_polar    % Drag coefficient
            % Third : alpha, cl, cd polar
            %processed_polars.al_alclcd_polar  % Angle of attack (deg)
            %processed_polars.cl_alclcd_polar  % Lift coefficient
            %processed_polars.cd_alclcd_polar  % Drag coefficent
            % Fourth: alpha, cm polar
            %processed_polars.al_alcm_polar    % Angle of attack (deg)
            %processed_polars.cm_alcm_polar    % Picth coefficient
            
            % % Take care of Cl polar
            if not(isempty(alcl_polar_dataset))
                % Extract data from cl_polar_dataset
                al_alcl_polar_raw       = alcl_polar_dataset(:,1);
                cl_alcl_polar_raw       = alcl_polar_dataset(:,2);
                
                % Sort cl_polar dataset
                [~, sort_index_cl_polar]=  sort(al_alcl_polar_raw);
                al_alcl_polar_sorted    = al_alcl_polar_raw(sort_index_cl_polar);
                cl_alcl_polar_sorted    = cl_alcl_polar_raw(sort_index_cl_polar);
                
                % Now find monotonous streak
                % Compute forward difference of Cl vs alpha
                dcl_alcl_polar_sorted   = diff(cl_alcl_polar_sorted);
                % Find indices of points where polar is decreasing
                index_dcl_decreasing    = find(dcl_alcl_polar_sorted < 0);
                % Find a point near zero angle of attack
                [~, index_al_min_alcl_polar_sorted] = min(al_alcl_polar_sorted.^2);
                % Find highest aoa below near-zero piont where polar is decreasing
                index_start_monotonous  = max(index_dcl_decreasing(index_dcl_decreasing < index_al_min_alcl_polar_sorted)) + 1;
                if isempty(index_start_monotonous)
                    index_start_monotonous = 1;
                end
                % Find lowest aoa above near-zero piont where polar is decreasing
                index_end_monotonous    = min(index_dcl_decreasing(index_dcl_decreasing > index_al_min_alcl_polar_sorted));
                if isempty(index_end_monotonous)
                    index_end_monotonous = length(al_alcl_polar_sorted);
                end
                
                % Now make monotonous streaks of cl polar
                al_alcl_polar_monotonous= al_alcl_polar_sorted(index_start_monotonous:index_end_monotonous);
                cl_alcl_polar_monotonous= cl_alcl_polar_sorted(index_start_monotonous:index_end_monotonous);
                
                % Now interpolants on monotonous streaks of cl polar
                al_alcl_polar_fun       = @(cl) interp1(cl_alcl_polar_monotonous , al_alcl_polar_monotonous, cl );
            else
                al_alcl_polar_sorted    = [];
                cl_alcl_polar_sorted    = [];
            end
            
            % % Take care of Cd polar
            if not(isempty(clcd_polar_dataset))
                % Now extract data from cd_polar_dataset (cl, cd)
                cl_clcd_polar_raw       = clcd_polar_dataset(:,1);
                cd_clcd_polar_raw       = clcd_polar_dataset(:,2);
                
                % Now sort cl, cd polar (not used for reinterpolation)
                [~, sort_index_clcd_polar] =  sort(cl_clcd_polar_raw);
                cl_clcd_polar_sorted    = cl_clcd_polar_raw(sort_index_clcd_polar);
                cd_clcd_polar_sorted    = cd_clcd_polar_raw(sort_index_clcd_polar);
            else
                cl_clcd_polar_sorted    = [];
                cd_clcd_polar_sorted    = [];
            end
            
            % % Reconstruct aoa of Cd polar with Cl polar
            if and(not(isempty(clcd_polar_dataset)), not(isempty(alcl_polar_dataset)))
                % Now reinterpolate alpha from monotonous streak of cl polar to find
                % alpha's of cd polar
                al_reinterpolated       = al_alcl_polar_fun(cl_clcd_polar_raw);
                % Filter NaNs out
                al_filtered             = al_reinterpolated(not(isnan(al_reinterpolated)));
                cl_filtered             = cl_clcd_polar_raw(not(isnan(al_reinterpolated)));
                cd_filtered             = cd_clcd_polar_raw(not(isnan(al_reinterpolated)));
                % Sort
                [~, index_al_filtered]  = sort(al_filtered);
                
                % Make alpha, cl, cd polar data
                al_alclcd_polar         = al_filtered(index_al_filtered);
                cl_alclcd_polar         = cl_filtered(index_al_filtered);
                cd_alclcd_polar         = cd_filtered(index_al_filtered);
            else
                al_alclcd_polar         = [];
                cl_alclcd_polar         = [];
                cd_alclcd_polar         = [];
            end
            
            % % Finally, take care of Cm polar
            if not(isempty(alcm_polar_dataset))
                % Extract data from cm_polar_dataset
                al_alcm_polar_raw       = alcm_polar_dataset(:,1);
                cm_alcm_polar_raw       = alcm_polar_dataset(:,2);
                
                % Sort cl_polar dataset
                [~, sort_index_cl_polar]=  sort(al_alcm_polar_raw);
                al_alcm_polar_sorted    = al_alcm_polar_raw(sort_index_cl_polar);
                cm_alcm_polar_sorted    = cm_alcm_polar_raw(sort_index_cl_polar);
            else
                al_alcm_polar_sorted    = [];
                cm_alcm_polar_sorted    = [];
            end
            
            % Now store data into a structure (using a template!)
            processed_polars = experimental_processed_polars();
            % First : alpha, cl polar
            processed_polars.al_alcl_polar    = al_alcl_polar_sorted;
            processed_polars.cl_alcl_polar    = cl_alcl_polar_sorted;
            % Second: cl, cd polar
            processed_polars.cl_clcd_polar    = cl_clcd_polar_sorted;
            processed_polars.cd_clcd_polar    = cd_clcd_polar_sorted;
            % Third : alpha, cl, cd polar
            processed_polars.al_alclcd_polar  = al_alclcd_polar;
            processed_polars.cl_alclcd_polar  = cl_alclcd_polar;
            processed_polars.cd_alclcd_polar  = cd_alclcd_polar;
            % Fourth: alpha, cm polar
            processed_polars.al_alcm_polar    = al_alcm_polar_sorted;
            processed_polars.cm_alcm_polar    = cm_alcm_polar_sorted;
        end
        
        function [processed_polars_re3e6clean, processed_polars_re6e6clean, processed_polars_re9e6clean] = process_polars_from_dataset(datasets)
            % % % Process Polars from Datasets
            % % Process polars for Re=3e6 datasets
            lg.msg(TR824_reader, 'Processing polars for Re=3e6');
            % Let us now extract a tripplet of datasets
            alcl_polar_dataset_re3e6clean     = datasets.A01_data;
            alcm_polar_dataset_re3e6clean     = datasets.A07_data;
            clcd_polar_dataset_re3e6clean     = datasets.B01_data;
            % And get some processed polars out of it
            [processed_polars_re3e6clean]     = TR824_reader.process_polar_triplet(alcl_polar_dataset_re3e6clean, alcm_polar_dataset_re3e6clean, clcd_polar_dataset_re3e6clean);
            % Add some metata to the processed polar
            processed_polars_re3e6clean.datataset_filnename = datasets.filename;
            processed_polars_re3e6clean.source              = datasets.source;
            
            % % Process polars for Re=6e6 datasets
            lg.msg(TR824_reader, 'Processing polars for Re=6e6');
            % Let us now extract a tripplet of datasets
            alcl_polar_dataset_re6e6clean     = datasets.A02_data;
            alcm_polar_dataset_re6e6clean     = datasets.A08_data;
            clcd_polar_dataset_re6e6clean     = datasets.B02_data;
            % And get some processed polars out of it
            [processed_polars_re6e6clean]     = TR824_reader.process_polar_triplet(alcl_polar_dataset_re6e6clean, alcm_polar_dataset_re6e6clean, clcd_polar_dataset_re6e6clean);
            % Add some metata to the processed polar
            processed_polars_re6e6clean.datataset_filnename = datasets.filename;
            processed_polars_re6e6clean.source              = datasets.source;
            
            % % Process polars for Re=9e6 datasets
            lg.msg(TR824_reader, 'Processing polars for Re=9e6');
            % Let us now extract a tripplet of datasets
            alcl_polar_dataset_re9e6clean     = datasets.A03_data;
            alcm_polar_dataset_re9e6clean     = datasets.A09_data;
            clcd_polar_dataset_re9e6clean     = datasets.B03_data;
            % And get some processed polars out of it
            [processed_polars_re9e6clean]     = TR824_reader.process_polar_triplet(alcl_polar_dataset_re9e6clean, alcm_polar_dataset_re9e6clean, clcd_polar_dataset_re9e6clean);
            % Add some metata to the processed polar
            processed_polars_re9e6clean.datataset_filnename = datasets.filename;
            processed_polars_re9e6clean.source              = datasets.source;
        end
        
        function [polar_conditions_re3e6clean, polar_conditions_re6e6clean, polar_conditions_re9e6clean] = make_polar_conditions()
            % % % Make experimental polar conditions
            % % Make experimental polar conditions for Re=3e6 clean polars
            lg.msg(TR824_reader, 'Making polar conditions for Re=3e6');
            % Create structure (object) from template
            polar_conditions_re3e6clean = experimental_polar_conditions();
            % Fill it in
            polar_conditions_re3e6clean.M        = 0;     % Mach number
            polar_conditions_re3e6clean.Re       = 3e6;   % Reynolds Number
            polar_conditions_re3e6clean.N_crit   = 9;     % Critical Amplification Factor of Tollmien-Schlichted waves (Van Ingen e^N method)
            polar_conditions_re3e6clean.xtr_top  = 0.99;  % x/c of forced transition point for top side (1 = free)
            polar_conditions_re3e6clean.xtr_bot  = 0.99;  % x/c of forced transition point for bottom side (1 = free)
            
            % % Make experimental polar conditions for Re=6e6 clean polars
            lg.msg(TR824_reader, 'Making polar conditions for Re=6e6');
            % Create structure (object) from template
            polar_conditions_re6e6clean = experimental_polar_conditions();
            % Fill it in
            polar_conditions_re6e6clean.M        = 0;     % Mach number
            polar_conditions_re6e6clean.Re       = 6e6;   % Reynolds Number
            polar_conditions_re6e6clean.N_crit   = 9;     % Critical Amplification Factor of Tollmien-Schlichted waves (Van Ingen e^N method)
            polar_conditions_re6e6clean.xtr_top  = 0.99;  % x/c of forced transition point for top side (1 = free)
            polar_conditions_re6e6clean.xtr_bot  = 0.99;  % x/c of forced transition point for bottom side (1 = free)
            
            % % Make experimental polar conditions for Re=9e6 clean polars
            lg.msg(TR824_reader, 'Making polar conditions for Re=9e6');
            % Create structure (object) from template
            polar_conditions_re9e6clean = experimental_polar_conditions();
            % Fill it in
            polar_conditions_re9e6clean.M        = 0;     % Mach number
            polar_conditions_re9e6clean.Re       = 9e6;   % Reynolds Number
            polar_conditions_re9e6clean.N_crit   = 9;     % Critical Amplification Factor of Tollmien-Schlichted waves (Van Ingen e^N method)
            polar_conditions_re9e6clean.xtr_top  = 0.99;  % x/c of forced transition point for top side (1 = free)
            polar_conditions_re9e6clean.xtr_bot  = 0.99;  % x/c of forced transition point for bottom side (1 = free)
        end
        
        function [result_flag_re3e6clean, result_flag_re6e6clean, result_flag_re9e6clean, ...
                  EC_re3e6clean         , EC_re6e6clean         , EC_re9e6clean         ] = make_and_add_experimental_cases(EDB, airfoil_description, polar_conditions_re3e6clean, processed_polars_re3e6clean, polar_conditions_re6e6clean, processed_polars_re6e6clean, polar_conditions_re9e6clean, processed_polars_re9e6clean)
            % Make experimental cases (and add them to database)
            % For Re = 3e6
            EC_re3e6clean = experimental_case(EDB);
            EC_re3e6clean.load_airfoil_description(airfoil_description)
            EC_re3e6clean.load_polar_conditions(polar_conditions_re3e6clean);
            EC_re3e6clean.load_processed_polars(processed_polars_re3e6clean);
            EC_re3e6clean.standardize_polars();
            result_flag_re3e6clean = add_experimental_case(EDB, EC_re3e6clean);
            
            % For Re = 6e6
            EC_re6e6clean = experimental_case(EDB);
            EC_re6e6clean.load_airfoil_description(airfoil_description)
            EC_re6e6clean.load_polar_conditions(polar_conditions_re6e6clean);
            EC_re6e6clean.load_processed_polars(processed_polars_re6e6clean);
            EC_re6e6clean.standardize_polars();
            result_flag_re6e6clean = add_experimental_case(EDB, EC_re6e6clean);
            
            % For Re = 9e6
            EC_re9e6clean = experimental_case(EDB);
            EC_re9e6clean.load_airfoil_description(airfoil_description)
            EC_re9e6clean.load_polar_conditions(polar_conditions_re9e6clean);
            EC_re9e6clean.load_processed_polars(processed_polars_re9e6clean);
            EC_re9e6clean.standardize_polars();
            result_flag_re9e6clean = add_experimental_case(EDB, EC_re9e6clean);
            
        end
        
        function [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = load_make_and_add_experimental_cases(EDB, airfoil_name, airfoil_filename, dataset_filename)
            % Load Make and Add Experimental Cases to EDB database
            
            % Load airfoil Description
            airfoil_description = EDB.SF.load_airfoil_description_from_file(airfoil_filename, airfoil_name);
            
            % Load datasets
            [datasets] = TR824_reader.read_datasets_from_file(dataset_filename);
            
            % Process Polars from Datasets
            [processed_polars_re3e6clean, processed_polars_re6e6clean, processed_polars_re9e6clean] = TR824_reader.process_polars_from_dataset(datasets);
            
            % Make experimental polar conditions
            [polar_conditions_re3e6clean, polar_conditions_re6e6clean, polar_conditions_re9e6clean] = TR824_reader.make_polar_conditions();
            
            % Make experimental cases (and add them to database)
            [result_flag_re3e6clean, result_flag_re6e6clean, result_flag_re9e6clean, ...
                EC_re3e6clean         , EC_re6e6clean         , EC_re9e6clean ] = ...
                TR824_reader.make_and_add_experimental_cases(EDB, airfoil_description, ...
                polar_conditions_re3e6clean, processed_polars_re3e6clean, ...
                polar_conditions_re6e6clean, processed_polars_re6e6clean, ...
                polar_conditions_re9e6clean, processed_polars_re9e6clean); %#ok<ASGLU>
        end
    end
end

