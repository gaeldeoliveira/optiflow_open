function sc = cost_function_RMSdCpdX_at_CL(parameters , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors

%% Input extraction
% Extract parameters from parameters structure
cl_i = parameters.cl_i;
w_i  = parameters.w_i;

% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{1};
ap = experiment_result.ap;

try
    if strcmp(class(ap) , 'aerodynamic_polar') %#ok<STISA>
        % If polar is valid
        if ~isempty(ap.bldata)
            % If bldata was properly read
            bldata = ap.bldata;
            
            %% Processing of Bldata
            xun_top_range   = bldata.xun_top_range;
            xun_bot_range   = bldata.xun_bot_range;
            cpx_top_rfun    = @(aoa)      bldata.cpx_fun_top(aoa, xun_top_range);
            cpx_bot_rfun    = @(aoa)      bldata.cpx_fun_bot(aoa, xun_bot_range);
            dcpx_dx_top_rfun= @(aoa)      gradient(cpx_top_rfun(aoa), xun_top_range);
            dcpx_dx_bot_rfun= @(aoa)      gradient(cpx_bot_rfun(aoa), xun_bot_range);
            %dcpx_dx_top_fun = @(aoa, xun) interp1(xun_top_range, dcpx_dx_top_rfun(aoa), xun);
            %dcpx_dx_bot_fun = @(aoa, xun) interp1(xun_bot_range, dcpx_dx_bot_rfun(aoa), xun);
            
            RMS_dcpx_dx_fun = @(aoa)      sqrt(trapz(xun_top_range, dcpx_dx_top_rfun(aoa).^2)) + sqrt(trapz(xun_bot_range, dcpx_dx_bot_rfun(aoa).^2));
            
            
            %% Processing
            
            % Allocate Variables
            alpha_i = zeros(size(cl_i));
            dcp_i   = zeros(size(cl_i));
            
            % Get cl range bounds
            ub_cl_range = max(ap.cl_vector_cl_range);
            lb_cl_range = min(ap.cl_vector_cl_range);
            width_cl_range = ub_cl_range - lb_cl_range;
            % And alpha of cl max and min
            alpha_cl_max = ap.alpha_cl(max(ap.cl_vector_cl_range));
            alpha_cl_min = ap.alpha_cl(min(ap.cl_vector_cl_range));
            
            
            % Find AOA's corresponding to desired Cl's
            for n_cl = 1:length(cl_i)
                alpha_i(n_cl) = ap.alpha_cl(cl_i(n_cl));
                
                if ~isnan(alpha_i)
                    % Cl range includes desired angle of attack
                    
                    % Find dCp's
                    dcp_i(n_cl) = RMS_dcpx_dx_fun(alpha_i(n_cl));
                else
                    % If Cl range does not cover desired cost function then return a
                    % valid value based on a penalty function
                    
                    if cl_i(n_cl) > ub_cl_range
                        % The requested cl is above the maximum one
                        
                        % Make penalty factor based on relative distance of
                        % requested Cl to achieved Cls of airfoil
                        relative_offset = (cl_i(n_cl)-ub_cl_range) / width_cl_range;
                        penalty_factor = 1 + abs(relative_offset);
                        
                        % Find dCp at maximum Cl
                        dcp_cl_max = RMS_dcpx_dx_fun(alpha_cl_max);
                        % Assign this value multiplied by penalty factor
                        dcp_i(n_cl) = dcp_cl_max * penalty_factor;
                        
                        disp(['Extrapolated for Higher Cls than those in range' , num2str(cl_i(n_cl))]);
                    end
                    
                    if cl_i(n_cl) < lb_cl_range
                        % If the requested cl is below the minimum one, just make a
                        % penalty factor and use the cl_cd at the lowest cl in range
                        
                        % Determine penalty factor
                        %   Follow same philosophy as for cl above range
                        relative_offset = (lb_cl_range-cl_i(n_cl)) / width_cl_range;
                        penalty_factor = 1 + abs(relative_offset);
                        
                        % Find dCp at minimum Cl
                        dcp_cl_min = RMS_dcpx_dx_fun(alpha_cl_min);
                        % Assign this value multiplied by penalty factor
                        dcp_i(n_cl) = dcp_cl_min * penalty_factor;
                        
                        disp(['Extrapolated for Lower Cls than those in range' , num2str(cl_i(n_cl))]);
                    end
                end
                
                % Weight Sum dCpdX values at different Cl's
                sc = sum(dcp_i  .* w_i);
                                
            end
        else
            % If bldata was unavailable
            sc = Inf;
        end        
        % Rule out "ultra-high performance" unphysical solutions (L/D > 300)
        if max(abs(ap.raw_data.cl ./ ap.raw_data.cd)) > 300
            sc = Inf;
        end        
    else
        % If polar is invalid!
        sc = Inf;
    end
    
    
    % Make check to verify that sides do not cross
    valid_shape = valid_shape_test(experiment_result);
    if not(valid_shape)
        sc = Inf;
    end
catch ME
    % Something went wrong! Display it:
    disp(['While evaluating: ' num2str(experiment_result.x)])
    disp(['Error occured in: ' mfilename()])
    disp(['Error identifier: ' ME.identifier])
    disp(['Error call stack:']);
    for is=1:length(ME.stack); disp(ME.stack(is)); end;
    % And assign invalidity value
    sc = Inf;    
end