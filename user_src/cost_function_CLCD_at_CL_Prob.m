function sc = cost_function_CLCD_at_CL(parameters , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors

%% Input extraction
% Extract parameters from parameters structure
cl_i0 = parameters.cl_i0;       % Reference Cl

% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{1};
ap = experiment_result.ap;

%% Extract last coordinate (equal to turbulence intensity)
x = experiment_result.x;
sigma_TI = x(end);
%
cl_i     = cl_i0 + linspace(-1,1,10) * sigma_TI   ;
w_i      =         ones(size(cl_i)) / length(cl_i);


%% Now run usual clcd at cl cost function


if strcmp(class(ap) , 'aerodynamic_polar')
    % If polar is valid
    
    %% Processing
    
    % Allocate Variables
    alpha_i = zeros(size(cl_i));
    cd_i    = zeros(size(cl_i));
    cl_cd   = zeros(size(cl_i));
    
    % Get cl range bounds
    ub_cl_range = max(ap.cl_vector_cl_range);
    lb_cl_range = min(ap.cl_vector_cl_range);
    width_cl_range = ub_cl_range - lb_cl_range;
    
    % Find AOA's corresponding to desired Cl's
    for n_cl = 1:length(cl_i)
        alpha_i(n_cl) = ap.alpha_cl(cl_i(n_cl));
        
        if ~isnan(alpha_i)
            % Cl range includes desired angle of attack
            
            % Find Cd's
            cd_i(n_cl) = ap.cd_alpha(alpha_i(n_cl));
            
            % Get Cl/Cds
            cl_cd(n_cl) = cl_i(n_cl) ./ cd_i(n_cl);
        else
            % If Cl range does not cover desired cost function then return a
            % valid value based on a penalty function
            
            if cl_i(n_cl) > ub_cl_range
                % The requested cl is above the maximum one
                
                % Determine an equivalent angle of attack offset based on
                % cl = 2*pi*alpha
                % alpha = cl / (2*pi())       <rad>
                % alpha = 90 / (pi()^2) * cl  <deg>
                
                cl_offset_absolute = (cl_i(n_cl) - ub_cl_range);
                alpha_offset_equivalent = 90/pi()^2 * cl_offset_absolute;
                
                % Damp out this formula to account for probable stall by
                % dubling alpha offset
                
                alpha_offset_equivalent = alpha_offset_equivalent * 2;
                
                %Determine Equivalent Angle of Attack
                alpha_cl_ub =  ap.alpha_cl(ub_cl_range);
                alpha_equiv = alpha_cl_ub + alpha_offset_equivalent;
                
                % And Cl and Cd at equivalent angle of attack
                cl_alpha_equiv = ap.cl_alpha(alpha_equiv);
                cd_alpha_equiv = ap.cd_alpha(alpha_equiv);
                
                % Determine cl_cd at equivalent angle of attack
                cl_cd(n_cl) = cl_alpha_equiv / cd_alpha_equiv;
                
                % Now impose a penalty factor to the whole thing for violating
                % constraints
                relative_offset = (cl_i(n_cl) - ub_cl_range) / width_cl_range;
                penalty_factor = 1 + relative_offset;
                
                % Apply Penalty Factor
                cl_cd(n_cl) = cl_cd(n_cl) / penalty_factor;
                
                disp(['Extrapolated for Higher Cls than those in range' , num2str(cl_i(n_cl))]);
            end
            
            if cl_i(n_cl) < lb_cl_range
                % If the requested cl is below the minimum one, just make a
                % penalty factor and use the cl_cd at the lowest cl in range
                
                % Determine alpha, cd and cl/cd at lowest cl in range
                alpha_cl_lb = ap.alpha_cl(lb_cl_range);
                cd_i(n_cl) = ap.cd_alpha(alpha_cl_lb);
                cl_cd(n_cl) = lb_cl_range / cd_i(n_cl);
                
                % Determine penalty factor
                %   Follow same philosophy as for cl above range
                relative_offset = (lb_cl_range-cl_i(n_cl)) / width_cl_range;
                penalty_factor = 1 + relative_offset;
                
                % Apply Penalty Factor
                cl_cd(n_cl) = cl_cd(n_cl) / penalty_factor;
                disp(['Extrapolated for Lower Cls than those in range' , num2str(cl_i(n_cl))]);
            end
        end
        
        % Weight Sum Cl/Cd ratio's
        sc = sum(cl_cd  .* w_i);
        
        %% Old Version
        %     % Find AOA's corresponding to desired Cl's
        %     alpha_i = ap.alpha_cl(cl_i);
        %
        %     if ~isnan(alpha_i)
        %         % Find Cd's
        %         cd_i = ap.cd_alpha(alpha_i);
        %
        %         % Get Cl/Cds
        %         cl_cd = cl_i ./ cd_i;
        %         % Weight Sum Cl/Cd ratio's
        %         sc = sum(cl_cd  .* w_i);
        %     else
        %         % If Cl range does not cover desired cost function then return NaN
        %         sc = NaN;
        %     end
        
    end
    
    % Rule out "ultra-high performance" unphysical solutions (L/D > 300)
    if max(abs(ap.raw_data.cl ./ ap.raw_data.cd)) > 300
        sc = 0;
    end
else
    % If polar is invalid!
    sc = 0;
end

% Make check to verify that sides do not cross
valid_shape = valid_shape_test(experiment_result);
if not(valid_shape)
   sc = 0; 
end

% Make additional check to verify that L/D did not enter into unreasonable
% branch