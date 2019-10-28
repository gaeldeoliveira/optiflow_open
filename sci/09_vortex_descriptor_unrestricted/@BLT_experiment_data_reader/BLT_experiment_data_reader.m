classdef BLT_experiment_data_reader < handle
    %BLT_experiment_data_reader is a handle class that reads, stores and
    % interpolates data from a the Low Speed Lab (LSL) Boundary Layer
    % Tunnel (BLT) Particle Image Velocimetry (PIV) experiments on a flat
    % plate w/wo VGs conducted by Baldacchino and Ragni
    
    properties
        data_folder = 'piv_data/';  % Default folder containing PIV data files
        fname_root                  % Root of PIV data file name (e.g. ZPG_yaw0_VGrect_d or ZPG_yaw0_base_d)
        fname_ext   = '.dat';       % Default extension of PIV data files
        
        d_list_flag                 % [bool] array of flags for filtering streamwise stances (1=filter, 0=do not filter at this x_over_h, as given in d_list)
        z_min_list                  % [m   ] array of filtering windows corner coordinates (left bound) (values only matter for points where d_list_flag is non-zero!)
        z_max_list                  % [m   ] array of filtering windows corner coordinates (right bound) (values only matter for points where d_list_flag is non-zero!)
        y_min_list                  % [m   ] array of filtering windows corner coordinates (bottom bound) (values only matter for points where d_list_flag is non-zero!)
        y_max_list                  % [m   ] array of filtering windows corner coordinates (top bound) (values only matter for points where d_list_flag is non-zero!)
        
        h                           % [m   ] Height of vane, used as scaling for x-direction
        S                           % [m   ] Half-spacing between VG pairs
        delta_ref                   % [m   ] Reference boundary layer height (needs not match BL height exactly)
        
        d_list                      % Line array of distances from TE in h units;
        data_struct_cell            % Cell array of piv data structures
    end
    
    methods
        function ED = BLT_experiment_data_reader(fname_root, d_list, h)
            % Intronize Inputed fields
            ED.fname_root = fname_root;
            ED.d_list     = d_list;
            ED.h          = h;
            
            ED.read_data();
        end
        
        function set_filters(ED, d_list_flag, z_min_list, z_max_list, y_min_list, y_max_list)
            % Intronize Inputed fields
            ED.d_list_flag = d_list_flag;
            ED.z_min_list  = z_min_list;
            ED.z_max_list  = z_max_list;
            ED.y_min_list  = y_min_list;
            ED.y_max_list  = y_max_list;
        end
        
        function read_data(ED)
            % Create Cell array of data_struct structures containing PIV
            % data for each downstream stance
            ED.data_struct_cell = cell(size(ED.d_list));
            
            % Read throught all steps
            for nd_list = 1:length(ED.d_list)
                % Get current distance identifier
                nd = ED.d_list(nd_list);
                % Make Filename
                if nd > 0
                    filename    = [ED.data_folder  ED.fname_root         , num2str(nd), ED.fname_ext ];
                else
                    filename    = [ED.data_folder  ED.fname_root(1:end-1), 'o'        , ED.fname_ext ];
                end
                ED.data_struct_cell{nd_list} = read_dat_file_into_data_struct(filename);
            end
        end
        
        function filter_data(ED)
            for nd_list = 1:length(ED.d_list)
                if not(ED.d_list_flag(nd_list) == false)
                    z_min = ED.z_min_list(nd_list);
                    z_max = ED.z_max_list(nd_list);
                    y_min = ED.y_min_list(nd_list);
                    y_max = ED.y_max_list(nd_list);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'u'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'v'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'w'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'uvw_length', z_min, z_max, y_min, y_max);
                    
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'u_rms'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'v_rms'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'w_rms'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                      'uvw_rms_length'  , z_min, z_max, y_min, y_max);
                    
                  
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'rst_zy'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'rst_zx'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'rst_yx'         , z_min, z_max, y_min, y_max);
                    
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'rst_zz'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'rst_yy'         , z_min, z_max, y_min, y_max);
                    ED.data_struct_cell{nd_list} = filter_data_struct_field( ED.data_struct_cell{nd_list}, ...
                        'rst_xx'         , z_min, z_max, y_min, y_max); 
                    
                    % Update interpolants
                    ED.data_struct_cell{nd_list}.u_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.u, z , y);
                    ED.data_struct_cell{nd_list}.v_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.v, z , y);
                    ED.data_struct_cell{nd_list}.w_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.w, z , y);
                    ED.data_struct_cell{nd_list}.interp_fun     = @(field,z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, field, z , y);
                end
            end
        end
            
        function set_crossflow_scaling_variables(ED, S, delta_ref)
            ED.S = S;
            ED.delta_ref = delta_ref;
        end
        
        function [u, v, w] = interp_velocities(ED, x_over_h, y_over_delta_ref, z_over_S)
            % Extrap with constant to sides in x (by fudging x_over_h)
            x_over_h = max([x_over_h, min(ED.d_list)]);
            x_over_h = min([x_over_h, max(ED.d_list)]);
            
            % Reconstruct dimensional plane variables
            y = y_over_delta_ref * ED.delta_ref;
            z = z_over_S         * ED.S;
            
            % Find index of experimental data just downstream (above) and upstream (below)
            nd_list_above = min(find(ED.d_list - x_over_h >= 0)); %#ok<MXFND>
            nd_list_below = max(find(ED.d_list - x_over_h <= 0)); %#ok<MXFND>
            
            % Find distances of just downstream (above) and upstream (below) stances
            x_over_h_above = ED.d_list(nd_list_above);
            x_over_h_below = ED.d_list(nd_list_below);
            
            % Interpolate Fields in Plane, just downstream (above) and upstream (below)
            u_above = ED.data_struct_cell{nd_list_above}.u_fun(z,y);
            u_below = ED.data_struct_cell{nd_list_below}.u_fun(z,y);
            
            v_above = ED.data_struct_cell{nd_list_above}.v_fun(z,y);
            v_below = ED.data_struct_cell{nd_list_below}.v_fun(z,y);
            
            w_above = ED.data_struct_cell{nd_list_above}.w_fun(z,y);
            w_below = ED.data_struct_cell{nd_list_below}.w_fun(z,y);
            
            % Mingle in Streamwise direction
            if not(nd_list_above == nd_list_below)
                u       = u_below + (u_above - u_below) / (x_over_h_above - x_over_h_below) * (x_over_h-x_over_h_below);
                v       = v_below + (v_above - v_below) / (x_over_h_above - x_over_h_below) * (x_over_h-x_over_h_below);
                w       = w_below + (w_above - w_below) / (x_over_h_above - x_over_h_below) * (x_over_h-x_over_h_below);
            else
                u       = u_below;
                v       = v_below;
                w       = w_below; 
            end
        end
        
        function [z, y, u] = interp_u_mesh(ED, x_over_h)
           % Extrap with constant to sides in x (by fudging x_over_h)
            x_over_h = max([x_over_h, min(ED.d_list)]);
            x_over_h = min([x_over_h, max(ED.d_list)]);
            
            % Find index of experimental data just downstream (above) and upstream (below)
            nd_list_above = min(find(ED.d_list - x_over_h >= 0)); %#ok<MXFND>
            nd_list_below = max(find(ED.d_list - x_over_h <= 0)); %#ok<MXFND>
            
            % Find distances of just downstream (above) and upstream (below) stances
            x_over_h_above = ED.d_list(nd_list_above);
            x_over_h_below = ED.d_list(nd_list_below);
            
            % Interpolate Fields in Plane, just downstream (above) and upstream (below)
            z_above = ED.data_struct_cell{nd_list_above}.z;
            z_below = ED.data_struct_cell{nd_list_below}.z;
            
            y_above = ED.data_struct_cell{nd_list_above}.y;
            y_below = ED.data_struct_cell{nd_list_below}.y;
            
            u_above = ED.data_struct_cell{nd_list_above}.u;
            u_below = ED.data_struct_cell{nd_list_below}.u;
            
            % Mingle in Streamwise direction
            if not(nd_list_above == nd_list_below)
                try
                    % Works whenever grids have the same sizes in all
                    % experimental snapshots (which is the case for ZPG,
                    % but not for APG)
                    z       = z_below + (z_above - z_below) / (x_over_h_above - x_over_h_below) * (x_over_h-x_over_h_below);
                    y       = y_below + (y_above - y_below) / (x_over_h_above - x_over_h_below) * (x_over_h-x_over_h_below);
                    u       = u_below + (u_above - u_below) / (x_over_h_above - x_over_h_below) * (x_over_h-x_over_h_below);
                catch
                    % If subsequent grids do not have the same size, stick
                    % to the one from below (later on assess oportunities for
                    % reinterpolation!)
                    z       = z_below;
                    y       = y_below;
                    u       = u_below;
                end
            else
                z       = z_below;
                y       = y_below; 
                u       = u_below;
            end 
        end
        
        function [dstr, theta, Hk, target_y, u_target_y] = integral_parameters(ED, x_over_h, z_over_S, y_over_delta_ref_range, u_inf)
            % Thicknesses are returned in [m   ], that is, rescaled from delta_ref units!
            
            % Reinterpolate through data
            target_y = linspace(0,y_over_delta_ref_range, 400)';
            u_target_y  = ED.interp_velocities(x_over_h , target_y, z_over_S) / u_inf;
            
            % Construct Integrands
            dstr_integrand =             (1 - u_target_y );
            theta_integrand = u_target_y  .* (1 - u_target_y );
            
            % Integrate Thicknesses with Trapezoidal Rule
            dstr = trapz(target_y, dstr_integrand) * ED.delta_ref;
            theta = trapz(target_y, theta_integrand) * ED.delta_ref;
            
            % Make Shape Factor
            Hk    = dstr / theta;
        end
        
        function extrapolate_u_fields(ED, u_inf, nu_inf)
            % Define number of points on added profile
            N_added_stances = 5;
            
            % Select current frame (this will become a loop counter)
            for nd_list = 1:length(ED.d_list)
                
                % Get boundary layer parameters for this frame
                target_z = linspace(-1,1);          % List of z/S stances
                Hk_z     = zeros(size(target_z ));  % Shape factor at z/S stance
                theta_z  = zeros(size(target_z ));  % Momentum thicknes at z/S stance (over delta_ref)
                for n_z = 1:length(target_z)
                    [~, theta_z(n_z), Hk_z(n_z)] = ED.integral_parameters(ED.d_list(nd_list),  target_z(n_z), 2.0, u_inf);
                end
                rt_z     = theta_z*u_inf/nu_inf;    % Reynolds theta at z/S stance
                
                % Take reference values as averages
                rt_mean    = mean(rt_z);
                Hk_mean    = mean(Hk_z);
                theta_mean = mean(theta_z);
                
                % Make high fidelity swafford profile
                SP                      = swafford_profile();
                SP.update_hk_rt_pair(Hk_mean, rt_mean);
                target_y_hf             =  linspace(0,2,4000);
                u_target_y_hf_swafford  = SP.evaluate_profile(target_y_hf * ED.delta_ref/ theta_mean);
                
                % Now get original u and y _meshes
                y_mesh_original_scaled = ED.data_struct_cell{nd_list}.y / ED.delta_ref;
                u_mesh_original_scaled = ED.data_struct_cell{nd_list}.u / u_inf;
                
                % And allocate space for extended u mesh
                u_mesh_extrapolated    = zeros(size(u_mesh_original_scaled));
                % And for offset vector
                y_offset_extrapolation = zeros(1,size(y_mesh_original_scaled, 2));
                
                % Now make an extrapolation for each row (this will become a nested loop)
                for n_z_index = 1:size(y_mesh_original_scaled, 2)
                    
                    % Extract current row
                    y_line_original_scaled = y_mesh_original_scaled(:,n_z_index);
                    u_line_original_scaled = u_mesh_original_scaled(:,n_z_index);
                    
                    % Reorder for extrapolation
                    % [target_y_l , sort_index]  = sort(y_line_original_scaled);
                    % u_target_y_l               = u_line_original_scaled(sort_index);
                    
                    % Reorder for extrapolation (including zero filtering)
                    % Filter
                    y_line_original_scaled_filtered = y_line_original_scaled(u_line_original_scaled>eps());
                    u_line_original_scaled_filtered = u_line_original_scaled(u_line_original_scaled>eps());
                    % Then reorder
                    [target_y_l , sort_index]  = sort(y_line_original_scaled_filtered);
                    u_target_y_l               = u_line_original_scaled_filtered(sort_index);
                    
                    % Add struct cell
                    % Find starting velocity of profile
                    u_target_y_l_start = u_target_y_l(1);
                    % Equivalent corresponding height on Swafford profile
                    y_l_start_offset = interp1(u_target_y_hf_swafford, target_y_hf, u_target_y_l_start);
                    % Make vector of added stances in y
                    y_l_added_stances = y_l_start_offset * linspace(0,1,N_added_stances)';
                    % Make vector of added stances in u
                    u_l_added_stances = interp1(target_y_hf, u_target_y_hf_swafford, y_l_added_stances);
                    % Make extended vector of stances in y
                    y_l_extended      = [y_l_added_stances(1:end-1); target_y_l(:) + y_l_start_offset];
                    % Make extended vector of stances in u
                    u_l_extended      = [u_l_added_stances(1:end-1); u_target_y_l(:)];
%                     % Now compare results
%                     plot(u_l_extended, y_l_extended); hold on
%                     plot(u_target_y_l, target_y_l  ); grid on
                    % Now re-interpolate onto extrapolated u_mesh (unscaled)
                    u_mesh_extrapolated(:,n_z_index) = interp1(y_l_extended, u_l_extended, y_line_original_scaled) * u_inf;
                    
                    % And store offset onto reference vector
                    y_offset_extrapolation(n_z_index) = y_l_start_offset;
                end
                
                % Feed result back into extrapolated mesh!
                ED.data_struct_cell{nd_list}.u_mesh_extrapolated    = u_mesh_extrapolated;
                % And store offsets for other processing needs!
                ED.data_struct_cell{nd_list}.y_offset_extrapolation = y_offset_extrapolation * ED.delta_ref;
            end
            
        end
        
        function make_extrapolated_fields_default(ED)
            for nd_list = 1:length(ED.d_list)
                % Create backup
               ED.data_struct_cell{nd_list}.u_original = ED.data_struct_cell{nd_list}.u;
               % Copu extrapolated to default
               ED.data_struct_cell{nd_list}.u          = ED.data_struct_cell{nd_list}.u_mesh_extrapolated;
               
               % Update interpolants
               ED.data_struct_cell{nd_list}.u_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.u, z , y);
               ED.data_struct_cell{nd_list}.v_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.v, z , y);
               ED.data_struct_cell{nd_list}.w_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.w, z , y);
               ED.data_struct_cell{nd_list}.interp_fun     = @(field,z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, field, z , y);
            end
        end
        
        function make_original_fields_default(ED)
            for nd_list = 1:length(ED.d_list)
               % Copy original to default
               ED.data_struct_cell{nd_list}.u          = ED.data_struct_cell{nd_list}.u_original;
               
               % Update interpolants
               ED.data_struct_cell{nd_list}.u_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.u, z , y);
               ED.data_struct_cell{nd_list}.v_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.v, z , y);
               ED.data_struct_cell{nd_list}.w_fun          = @(z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, ED.data_struct_cell{nd_list}.w, z , y);
               ED.data_struct_cell{nd_list}.interp_fun     = @(field,z,y) interp2(ED.data_struct_cell{nd_list}.z, ED.data_struct_cell{nd_list}.y, field, z , y);
            end
        end
        
    end
    
end

