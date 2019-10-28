function sc = cost_function_CpFinal(parameters , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors

%% Input extraction
% Extract parameters from parameters structure

% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{1};
ap = experiment_result.ap;
pitch_angle = experiment_result.x(end);

if strcmp(class(ap) , 'aerodynamic_polar')
    % If polar is valid
    
    % find the appropriate thickness (of the NACA00XX airfoil) interpolated between the simulations 
    CldAlpha = (interp1(ap.raw_data.alpha , ap.raw_data.cl , 3, 'linear' , 0)...
                -interp1(ap.raw_data.alpha , ap.raw_data.cl , 0, 'linear' , 0))/3; 

    [MAXCL,IMAXCL]=max(ap.raw_data.cl);
    [MINCL,IMINCL]=min(ap.raw_data.cl);
    sc=(MAXCL-MINCL)./(ap.raw_data.alpha(IMAXCL(1))-ap.raw_data.alpha(IMINCL(end)));

    if sc<0.086
      sc = Inf;
    end        
    
else
    % If polar is invalid!
    sc = Inf;
end