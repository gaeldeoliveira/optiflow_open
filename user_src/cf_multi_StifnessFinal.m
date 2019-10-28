function sc = cf_multi_StifnessFinal( ~ , experiments_results)
% Calculate stiffness of airfoil
% Extract geometric properties from standard experiment input form
result_soiled  = experiments_results{1};
ap_soiled = result_soiled.ap;
result_clean   = experiments_results{2};
ap_clean = result_clean.ap;

alpha_range = 15;

coordinates = result_clean.coordinates;
COORD(:,1)=coordinates.tx;
COORD(:,2)=coordinates.tz;

AIRFOIL=calculate_structural_properties_airfoils(COORD,300);
sc=AIRFOIL.YBEND.Ixx_over_tsurf;

%%%% correct by trailing edge angle
% max_angle = 12;
% PENALTY=sqrt(max(AIRFOIL.ANGLE_TE,0)/max_angle);
% PENALTY=min(PENALTY,1);
% sc1=sc1*PENALTY;
% 
%    disp(['     STIFNESS = ' num2str(sc1) ])
%    disp(['     TE ANGLE = ' num2str(AIRFOIL.ANGLE_TE) ])
%    disp(['     STIFNESS_CORRECTED_TE = ' num2str(sc2) ])
% 
% sc=sc1;

%%%% thickness check + penalty
% if AIRFOIL.THICKNESS_MAX > 0.33
%     sc=sc*(0.33/(abs(AIRFOIL.THICKNESS_MAX)+0.1)).^100;
% end

%%% aerodynamic check
if and(strcmp(class(ap_clean) , 'aerodynamic_polar') , strcmp(class(ap_soiled) , 'aerodynamic_polar'))
    % If polar is valid
    [SCMAX1,ANGLE_SCMAX1,~]=calculate_slope_function(ap_clean,alpha_range,0);
    [~,~,SC_ANGLE_IN2]=calculate_slope_function(ap_soiled,alpha_range,ANGLE_SCMAX1);         
    
    SCMAX1=max(SCMAX1,0.0000001);
    SC_ANGLE_IN2=max(SC_ANGLE_IN2,0.0000001);
    
    SCMAX=1/(0.5/SCMAX1+0.5/SC_ANGLE_IN2);
else
    SCMAX = 1;
end

% aerodynamic penalty
if SCMAX < 7.5
    sc = sc*SCMAX/12;
end
   disp(['     STIFNESS CORRECTED TE and SLOPE = ' num2str(sc) ])
   disp(['     SLOPE VALUE = ' num2str(SCMAX) ])
 
end