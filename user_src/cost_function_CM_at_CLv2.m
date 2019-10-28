function sc = cost_function_CM_at_CLv2(parameters , experiments_results)
%function sc = sub_cost_function_CM_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted CM of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors
%   Fudging and penalty factors are taylored to handle 

%% Input extraction
% Extract parameters from parameters structure
cl_i = parameters.cl_i;
w_i  = parameters.w_i;

% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{1};
ap = experiment_result.ap;

if strcmp(class(ap) , 'aerodynamic_polar')
    % If polar is valid

%% Processing

% Allocate Variables
alpha_i = zeros(size(cl_i));
cm_i    = zeros(size(cl_i));
cm      = zeros(size(cl_i));

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
        cm_i(n_cl) = ap.cm_alpha(alpha_i(n_cl));
    
        % Get Cl/Cds
        cm(n_cl)   = cm_i(n_cl);
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
            % doubling alpha offset
            
            alpha_offset_equivalent = alpha_offset_equivalent * 2;
                        
            %Determine Equivalent Angle of Attack
            alpha_cl_ub =  ap.alpha_cl(ub_cl_range);
            alpha_equiv = alpha_cl_ub + alpha_offset_equivalent;
            
            % And Cl and Cd at equivalent angle of attack
            cl_alpha_equiv = ap.cl_alpha(alpha_equiv);
            cm_alpha_equiv = ap.cm_alpha(alpha_equiv);
    
            % Determine cl_cd at equivalent angle of attack
            cm(n_cl) = cl_alpha_equiv / cm_alpha_equiv;
            
            % Now impose a penalty factor to the whole thing for violating
            % constraints
            relative_offset = (cl_i(n_cl) - ub_cl_range) / width_cl_range;
            penalty_factor = 1 + relative_offset;           
            
            % Apply Penalty Factor
            cm(n_cl) = cm(n_cl) * penalty_factor;
            
            disp(['Extrapolated for Higher Cls than those in range' , num2str(cl_i(n_cl))]);
        end
        
        if cl_i(n_cl) < lb_cl_range
            % If the requested cl is below the minimum one, just make a
            % penalty factor and use the cl_cd at the lowest cl in range
            
            % Determine alpha, cd and cl/cd at lowest cl in range
            alpha_cl_lb = ap.alpha_cl(lb_cl_range);
            cm_i(n_cl)  = ap.cm_alpha(alpha_cl_lb);
            cm(n_cl)    = lb_cl_range / cm_i(n_cl);
            
            % Determine penalty factor
            %   Follow same philosophy as for cl above range
            relative_offset = (lb_cl_range-cl_i(n_cl)) / width_cl_range;
            penalty_factor = 1 + relative_offset;
            
            % Apply Penalty Factor
            cm(n_cl) = cm(n_cl) * penalty_factor;            
            disp(['Extrapolated for Lower Cls than those in range' , num2str(cl_i(n_cl))]);
        end
    end

    % Weight Sum Cl/Cd ratio's
    sc = sum(cm  .* w_i);
    
    % Fudge value to a maximum of unity (way above anything reasonable for
    % a Cm!) (this fudging is good when you want to limit CM, not
    % otherwise)
    if abs(sc) > 1
        sc = -1;
    end





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

else
    % If polar is invalid!
    sc = 0;
end