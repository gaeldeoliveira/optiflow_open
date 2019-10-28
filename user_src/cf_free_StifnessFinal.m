function sc = cf_StifnessFinal( ~ , experiments_results)
% Calculate stiffness of airfoil
% Extract geometric properties from standard experiment input form
exp_result   = experiments_results{1};
% parameters = experiment_result.x;
coordinates = exp_result.coordinates;
COORD(:,1)=coordinates.tx;
COORD(:,2)=coordinates.tz;

AIRFOIL=calculate_structural_properties_airfoils(COORD,300);
sc1=AIRFOIL.YBEND.Ixx_over_tsurf;

%%% correct by trailing edge angle
PENALTY=sqrt(max(AIRFOIL.ANGLE_TE,0)/8);
PENALTY=min(PENALTY,1);
sc2=sc1*PENALTY;

   disp(['     STIFNESS = ' num2str(sc1) ])
   disp(['     TE ANGLE = ' num2str(AIRFOIL.ANGLE_TE) ])
   disp(['     STIFNESS_CORRECTED_TE = ' num2str(sc2) ])

sc=sc2;
% thickness check + penalty
if AIRFOIL.THICKNESS_MAX > 0.33
    sc=sc2*(0.33/(abs(AIRFOIL.THICKNESS_MAX)+0.1)).^100;
end

% aerodynamic check
if strcmp(class(exp_result.ap) , 'aerodynamic_polar')
    [SCMAX,~,~] = calculate_slope_function(exp_result.ap,15,0);
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