% A very simple BEM code bogo2_alfa7 hook to GAM optimization
% Driver for driven script
%
% To run publishable results:
%   Increase BC.N_bins      to 40/60
%   Increase PopulationSize to  ~200
%   Increase MaxGenerations to  ~100 (strict minimum)
%
% Polar tensor issues
%   ../rotor_integration/airfoil_families/DU-WP/FFA_DUWP_S20E50_free.mat

% Add needed paths
addpath src

% % Design problem settings (common to study)
PB.make_optimization        = true  ;
PB.MaxGenerations           =  80   ;
PB.PopulationSize           = 240   ;
PB.ParetoFraction           = 0.4   ;
PB.UseParallel              =  true ;
PB.base_blade_geometry_file = '../rotor_integration/graph_digitization/inwind_moded_RWT.mat'  ;
PB.std_distribution_file    = '../rotor_integration/probabilistic/std_distribution_scaled.mat';
if PB.make_optimization == true
    PB.start_datestr            = ['DO_' , datestr(now,30)];
else
    PB.start_datestr            = ['DT_' , datestr(now,30)];
end

% % Make folder to store results
mkdir([PB.start_datestr])

% % Stat, FFA foils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/FFA/FFA_freeC.mat'     ;
PB.bool_prob_coef           = false;
PB.case_id                  = 'stat_FFA_free';
% Run and save
bogo2_alfa7_realistic_rotor_optimization_driven_v2

% % Stat, DUIW High glide airfoils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/DU-IW/high_glide_design_free.mat';
PB.bool_prob_coef           = false;
PB.case_id                  = 'stat_high_glide_free';
% Run and save
bogo2_alfa7_realistic_rotor_optimization_driven_v2

% % Stat, DUIW Smooth glide airfoils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/DU-IW/smooth_glide_design_free.mat';
PB.bool_prob_coef           = false;
PB.case_id                  = 'stat_smooth_glide_free';
% Run and save
bogo2_alfa7_realistic_rotor_optimization_driven_v2


% % Prob, FFA foils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/FFA/FFA_freeC.mat'     ;
PB.bool_prob_coef           = true;
PB.case_id                  = 'prob_FFA_free';
% Run and save
bogo2_alfa7_realistic_rotor_optimization_driven_v2

% % Prob, DUIW High glide airfoils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/DU-IW/high_glide_design_free.mat';
PB.bool_prob_coef           = true;
PB.case_id                  = 'prob_high_glide_free';
% Run and save
bogo2_alfa7_realistic_rotor_optimization_driven_v2

% % Prob, DUIW Smooth glide airfoils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/DU-IW/smooth_glide_design_free.mat';
PB.bool_prob_coef           = true;
PB.case_id                  = 'prob_smooth_glide_free';
% Run and save
bogo2_alfa7_realistic_rotor_optimization_driven_v2






