classdef shape_fit_cst < shape_fit
    %SHAPE_DEFINITION This class stores, aquires and converts airfoil shape
    %definitions
     
    properties
        
    end
    
    properties(SetAccess = private) 
        regenerated_coordinates
        z_te_upper
        z_te_lower
        
        x_list_fitted
        names_fitted
        consistency_report_fitted
        
%       parameters_down
%       parameters_suction
    end
    
    methods
        function SF = shape_fit_cst(shape_definition_handle, name)
            SF.shape_definition_handle = shape_definition_handle;
            SF.name = name;
            SF.id = now();
        end        
        function x = get_parameters_from_file(SF , fileToRead)
            disp(['Loading and fitting airfoil from file : ' fileToRead])
            SF.import_coordinates(fileToRead);
            SF.fit_parameters_to_coordinates();
            SF.update_parameter_vector();
            SF.regenerate_coordinates();
            x = SF.parameters.full_non_linear.full_export;
        end        
        function fit_parameters_to_coordinates(SF)
            % Call fitting function with class specific parameters
            % Requires adimensionalized coordinates to be defined and the 
            [SF.parameters.linear.upper , SF.z_te_upper , ...
                SF.parameters.linear.lower , SF.z_te_lower , ...
                SF.parameters.non_linear.upper , SF.parameters.non_linear.lower , ...
                SF.parameters.full_non_linear.upper , SF.parameters.full_non_linear.lower] ...
                = cst_fit_airfoil(SF.shape_definition_handle.parametrization_handles.upper , ...
                                  SF.shape_definition_handle.parametrization_handles.lower , ...
                                  SF.imported_coordinates.adimensionalized.tx_coordinates , ...
                                  SF.imported_coordinates.adimensionalized.tz_coordinates);                                          
        end     
        function x_list = obtain_airfoil_list_from_folder(SF , airfoil_subdir)
            
            % FIrst list folder contents with ls, and process them into
            % separate files
            airfoil_list_string_array=[];

            airfoil_list_string  = dir([airfoil_subdir '*.txt']);
 %           length(airfoil_list_string)
            for k1=1:length(airfoil_list_string)
             airfoil_list_string_array{k1} = airfoil_list_string(k1).name; %#ok<AGROW>
            end
            k2p=length(airfoil_list_string);
  
            airfoil_list_string  = dir([airfoil_subdir '*.dat']);
 %           length(airfoil_list_string)          
            for k1=1:(length(airfoil_list_string))
             airfoil_list_string_array{k1+k2p} = airfoil_list_string(k1).name; %#ok<AGROW>
            end

            k2p=length(airfoil_list_string);
  
            airfoil_list_string  = dir([airfoil_subdir '*.air']);
 %           length(airfoil_list_string)          
            for k1=1:(length(airfoil_list_string))
             airfoil_list_string_array{k1+k2p} = airfoil_list_string(k1).name; %#ok<AGROW>
            end
            
            
            
            % Now identify how many airfoils there are to load in the
            % folder
            s = size(airfoil_list_string_array);
            N_airfoils = s(2);
            
            % Preallocate memory accordingly
            x_list = cell(N_airfoils , 1);
            
            % clean formerly allocated stuff
            SF.names_fitted = [];
            SF.consistency_report_fitted = [];
            
            % Now fit all airfoils in folder
            for n_airfoil = 1:N_airfoils
                % Construct current file name and extract it from singleton cell array
                current_airfoil_file = [airfoil_subdir  airfoil_list_string_array{n_airfoil}];
   %             current_airfoil_file = current_airfoil_file{1};
                % Fit !
                x_list{n_airfoil} = SF.get_parameters_from_file(current_airfoil_file);
                SF.names_fitted{n_airfoil}  = airfoil_list_string_array(n_airfoil);
                if n_airfoil == 1
                    SF.consistency_report_fitted = SF.consistency_check;
                else
                    SF.consistency_report_fitted(n_airfoil)  = SF.consistency_check;
                end
            end
            % Done ! Just store into reference population array for later
            % use
            SF.x_list_fitted = cell2mat(x_list);
            
        end
        function regenerate_coordinates(SF)
            % Generates tz coordinates from parameters at same tx points as
            % where original airfoil is defined
            
            
            %% Separate Upper and Lower Sides
            
            tx_coordinates = SF.imported_coordinates.adimensionalized.tx_coordinates;
            % Find leading edge index            
            [~ , min_index] = min(tx_coordinates);

            % Generate Upper and Lower tx coordinate vectors
            tx_lower = tx_coordinates(min_index:end);
            tx_upper = tx_coordinates(1:min_index-1); % Reorder upper side to avoid confusions further down
            
            fieldlist = fieldnames(SF.parameters);
            
            for n_field = 1:length(fieldlist)
                
                % Extract Parameters from reference field
                parameters_upper = SF.parameters.(fieldlist{n_field}).upper;
                parameters_lower = SF.parameters.(fieldlist{n_field}).lower;
                
                % Now generate ordinate vectors
                tz_upper = SF.shape_definition_handle.tz_upper(tx_upper , parameters_upper , SF.z_te_upper);
                tz_lower = SF.shape_definition_handle.tz_lower(tx_lower , parameters_lower , SF.z_te_lower);
                
                % Concatenate Vectors
                tx = [tx_upper ; tx_lower];
                tz = [tz_upper ; tz_lower];
                
                % With each set of fitted parameters generate ordinates at
                % same abcissa as original (adimensionalized) data
                SF.regenerated_coordinates.(fieldlist{n_field}).tx = tx;
                SF.regenerated_coordinates.(fieldlist{n_field}).tz = tz;                
            end            
        end
        function update_parameter_vector(SF)                        
            fieldlist = fieldnames(SF.parameters);
            
            % Compose the full shape parameter vectors from, upper, lower
            % and te parameters
            for n_field = 1:length(fieldlist)
                SF.parameters.(fieldlist{n_field}).full = ...
                    [SF.parameters.(fieldlist{n_field}).upper , ...
                    SF.parameters.(fieldlist{n_field}).lower , ...
                    SF.z_te_upper];                                    
            end
            
            % The above parameters are in the form of a pure CST shape
            % definition, and so have no dummy parameters.
            % Make a version of the results padded with non shape elements
            % for exportation to constraint module. Non shape elements are
            % set to zero.
            dummy_add = zeros(1 , SF.shape_definition_handle.N_dummy_parameters);
            for n_field = 1:length(fieldlist)
                SF.parameters.(fieldlist{n_field}).full_export = ...
                    [ SF.parameters.(fieldlist{n_field}).full , ...
                    dummy_add];
            end
            
        end
        function airfoil_description = load_airfoil_description_from_file(SF, fileToRead, airfoil_name)
            % Create airfoil description from template
            airfoil_description                    = experimental_airfoil_description();
            % Read airfoil from file!
            airfoil_description.x                  = SF.get_parameters_from_file(fileToRead);
            % Store additional metadata!
            airfoil_description.N_dummy_parameters = SF.shape_definition_handle.N_dummy_parameters;
            % Standardized Coordinates
            airfoil_description.tx_coordinates = SF.imported_coordinates.adimensionalized.tx_coordinates;
            airfoil_description.tz_coordinates = SF.imported_coordinates.adimensionalized.tz_coordinates;
            % Raw Coordinates from file
            airfoil_description.tx_coordinates_raw  = SF.imported_coordinates.raw.tx_coordinates;
            airfoil_description.tz_coordinates_raw  = SF.imported_coordinates.raw.tz_coordinates;
            % Raw to standardized transformation parameters
            airfoil_description.theta               = SF.imported_coordinates.processing_parameters.theta;
            airfoil_description.scale_factor        = SF.imported_coordinates.processing_parameters.scale_factor;
            airfoil_description.translation_factor  = SF.imported_coordinates.processing_parameters.translation_factor;
            % Airfoil Filename
            airfoil_description.airfoil_filename    = fileToRead;
            airfoil_description.airfoil_name        = airfoil_name;
        end
        function consistency_report = consistency_check(SF)

            % Consistency Check

            tx_ref = SF.imported_coordinates.adimensionalized.tx_coordinates;
            tz_ref = SF.imported_coordinates.adimensionalized.tz_coordinates;
            
            tx_reg = SF.regenerated_coordinates.full_non_linear.tx;
            tz_reg = SF.regenerated_coordinates.full_non_linear.tz;
            
            
            % Separate upper and lower sides and construct interpolants and derivatives
            
            % Find leading edge index
            [~ , min_index] = min(tx_reg);
            
            % Generate Upper and Lower tx coordinate vectors
            tx_reg_lower_side = tx_reg(min_index:end);
            tz_reg_lower_side = tz_reg(min_index:end);
            
            tx_reg_upper_side = flipud(tx_reg(1:min_index)); % Reorder upper side to avoid confusions further down
            tz_reg_upper_side = flipud(tz_reg(1:min_index)); % Reorder upper side to avoid confusions further down
            
            
            % Find leading edge index
            [~ , min_index] = min(tx_ref);
            
            % Generate Upper and Lower tx coordinate vectors
            tx_ref_lower_side = tx_ref(min_index:end); %#ok<NASGU>
            tz_ref_lower_side = tz_ref(min_index:end);
            
            tx_ref_upper_side = flipud(tx_ref(1:min_index)); %#ok<NASGU> % Reorder upper side to avoid confusions further down
            tz_ref_upper_side = flipud(tz_ref(1:min_index)); % Reorder upper side to avoid confusions further down
            
            
            % Build interpolants and derivative estimators
            % DEPRECATE due to change in Matlab R2013A+
                % pp_reg_upper = interp1(tx_reg_upper_side , tz_reg_upper_side , 'pchip' , 'pp');
            % Switch to pchip directly (stricly equivalent)
            pp_reg_upper = pchip(tx_reg_upper_side , tz_reg_upper_side);
            pp_reg_upper_d1 = fnder(pp_reg_upper);
    
            % DEPRECATE due to change in Matlab R2013A+
                % pp_reg_lower = interp1(tx_reg_lower_side , tz_reg_lower_side , 'pchip' , 'pp');
            % Switch to pchip directly (stricly equivalent)
            pp_reg_lower = pchip(tx_reg_lower_side , tz_reg_lower_side);
            pp_reg_lower_d1 = fnder(pp_reg_lower);
            
            % Find derivative values at relevant places
            reg_upper_d1_val = ppval(pp_reg_upper_d1 , tx_reg_upper_side);
            reg_lower_d1_val = ppval(pp_reg_lower_d1 , tx_reg_lower_side);
            
            
            % Compute height differences everywhere
            
            dv_upper = (tz_ref_upper_side - tz_reg_upper_side);
            dn_upper = dv_upper .* cos(atan(reg_upper_d1_val));
            
            dv_lower = (tz_ref_lower_side - tz_reg_lower_side);
            dn_lower = dv_lower .* cos(atan(reg_lower_d1_val));
            
            consistency_report.coordinates.tx_reg = tx_reg;
            consistency_report.coordinates.tz_reg = tz_reg;
            consistency_report.coordinates.tx_ref = tx_ref;
            consistency_report.coordinates.tz_ref = tz_ref;
            
            consistency_report.dv_upper  = dv_upper;
            consistency_report.dn_upper  = dn_upper;
            consistency_report.dv_lower  = dv_lower;
            consistency_report.dn_lower  = dn_lower;
            
            consistency_report.Npoints_definition = length(tx_ref); 
            
            consistency_report.RMS_upper = norm(dn_upper) ./ length(dn_upper);
            consistency_report.RMS_lower = norm(dn_lower) ./ length(dn_lower);
            consistency_report.RMS_err   = sqrt(norm(dn_upper)^2 + norm(dn_lower)^2) / ( length(dn_upper) + length(dn_lower));
            
            consistency_report.RMS_err   = sqrt(norm(dv_upper)^2 + norm(dv_lower)^2) / ( length(dv_upper) + length(dv_lower));
            
            consistency_report.max_err   = max( max(dn_upper)  , max(dn_lower) );
            consistency_report.max_verr  = max( max(dv_upper)  , max(dv_lower) );
        end
        %%% CONSISTENCY CHECK HAS TO BE REWRITTEN%%%
    end
    % Restons Amants - Maxime le Forestier
end

