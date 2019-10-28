function sc = cf_multi_StifnessFinal_box( ~ , experiments_results)
% Calculate stiffness of airfoil
% Extract geometric properties from standard experiment input form
exp_result   = experiments_results{2};
% parameters = experiment_result.x;
coordinates = exp_result.coordinates;
COORD(:,1)=coordinates.tx;
COORD(:,2)=coordinates.tz;

%%
box.front = 0.2;        % chordwise location of front spar of box
box.back = 0.7;         % chordwise location of back spar of box
box.web_t = 1;          % thickness of box web plates (spars) in multiples of skin thickness

AIRFOIL=calculate_structural_properties_airfoils_box(COORD,300,box);
Ixx_box = AIRFOIL.YBEND.Ixx_box_skin;% + AIRFOIL.YBEND.Ixx_box_webs;          % both webs and skin is Ixx/t

sc1 = Ixx_box;

%% correct by trailing edge angle
max_angle = 8;
PENALTY=sqrt(max(AIRFOIL.ANGLE_TE,0)/max_angle);
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