function sc = cf_aero_free(~ , experiments_results)
% Cost function for multisimulation case (forced and free transition)
% ....

% Results for forced(soiled) and free(clean) transition are loaded from results
result_clean   = experiments_results{1};
ap_clean = result_clean.ap;

alpha_range = 15;

if strcmp(class(ap_clean) , 'aerodynamic_polar')
    % If polar is valid
    [SCMAX1,~,~]=calculate_slope_function(ap_clean,alpha_range,0);
    %[SCMAX2,ANGLE_SCMAX2,SC_ANGLE_IN2]=calculate_slope_function(ap_soiled,alpha_range,ANGLE_SCMAX1);         
    
    sc=max(SCMAX1,0.0000001);
    
    disp('---------------------NEW AIRFOIL----------------------')
    disp(['Cost function dcl is ' num2str(sc)])
else
    % If polar is invalid!
    sc = -Inf;                  % minus because of post function
    disp('invalid polar')

end

end


