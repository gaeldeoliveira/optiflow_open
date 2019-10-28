function sc = cf_multi_AeroSlope(~ , experiments_results)
% Cost function for multisimulation case (forced and free transition)
% ....

% Results for forced(soiled) and free(clean) transition are loaded from results
result_soiled  = experiments_results{1};
ap_soiled = result_soiled.ap;
result_clean   = experiments_results{2};
ap_clean = result_clean.ap;

alpha_range = 15;

if and(strcmp(class(ap_clean) , 'aerodynamic_polar') , strcmp(class(ap_soiled) , 'aerodynamic_polar'))
    % If polar is valid
    [SCMAX1,ANGLE_SCMAX1,SC_ANGLE_IN1]=calculate_slope_function(ap_clean,alpha_range,0);
    [SCMAX2,ANGLE_SCMAX2,SC_ANGLE_IN2]=calculate_slope_function(ap_soiled,alpha_range,ANGLE_SCMAX1);         
    
    SCMAX1=max(SCMAX1,0.0000001);
    SC_ANGLE_IN2=max(SC_ANGLE_IN2,0.0000001);
    
    sc=1/(0.5/SCMAX1+0.5/SC_ANGLE_IN2);
    
    disp('---------------------NEW AIRFOIL----------------------')
    disp(['Cost function dcl is ' num2str(sc)])
    disp(['Cost function clean  ' num2str(SCMAX1) ])
    disp(['Cost function soiled ' num2str(SC_ANGLE_IN2)])
    disp(['Pitch angle is       ' num2str(ANGLE_SCMAX1) ])
else
    % If polar is invalid!
    sc = -Inf;                  % minus because of post function
    disp('invalid polar')

end

end


