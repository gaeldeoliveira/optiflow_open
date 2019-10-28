% Add needed paths
addpath src

% % Design problem settings (common to study)
PB.make_optimization        = false  ;
PB.MaxGenerations           =  80   ;
PB.PopulationSize           = 224   ;
PB.ParetoFraction           = 0.4   ;
PB.UseParallel              =  true ;
PB.base_blade_geometry_file = '../rotor_integration/graph_digitization/inwind_moded_RWT.mat';
PB.std_distribution_file    = '../rotor_integration/probabilistic/std_distribution_scaled.mat'     ;
if PB.make_optimization == true
    PB.start_datestr            = ['DO_' , datestr(now,30)];
else
    PB.start_datestr            = ['DT_' , datestr(now,30)];
end
% % Restart
% PB.start_datestr = 'DO_20180908T005423';

% % Make folder to store results
% mkdir([PB.start_datestr])

% % Stat, FFA foils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/FFA/FFA_freeC.mat'     ;
PB.bool_prob_coef           = true;
% Initialize interpolators
bogo2_alfa7_realistic_polar_comparison_driven
% Store interpolators
BC_stat_FFA_free = BC;

% % Stat, DUWP S00E35 foils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/DU-WP/FFA_DUWP_S00E35_free.mat';
PB.bool_prob_coef           = true;
% Initialize interpolators
bogo2_alfa7_realistic_polar_comparison_driven
% Store interpolators
BC_stat_S00E35_free = BC;

% % Stat, DUWP S20E35 foils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/DU-WP/FFA_DUWP_S20E35_free.mat';
PB.bool_prob_coef           = true;
% Initialize interpolators
bogo2_alfa7_realistic_polar_comparison_driven
% Store interpolators
BC_stat_S20E35_free = BC;

% % Stat, DUWP S20E35 foils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/DU-WP/FFA_DUWP_S40E35_trip.mat';
PB.bool_prob_coef           = true;
% Initialize interpolators
bogo2_alfa7_realistic_polar_comparison_driven
% Store interpolators
BC_stat_S40E35_free = BC;



% Now plot
tc_plot = 0.21; 
re_plot = 9e6 ;
mu_plot = 0.8 ;
alpha_range = 0:0.2:20;
figure(501)
subplot(221)
plot(alpha_range, BC_stat_S00E35_free.cl_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot));hold on
plot(alpha_range, BC_stat_S20E35_free.cl_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot));hold on
plot(alpha_range, BC_stat_S40E35_free.cl_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot));hold on
legend('S00E35', 'S20E35', 'S40E35'); title('Static Cl');

subplot(223)
plot(alpha_range, BC_stat_S00E35_free.cl_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot) ./ BC_stat_S00E35_free.cd_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot));hold on
plot(alpha_range, BC_stat_S20E35_free.cl_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot) ./ BC_stat_S20E35_free.cd_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot));hold on
plot(alpha_range, BC_stat_S40E35_free.cl_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot) ./ BC_stat_S40E35_free.cd_of_alpha_deg_and_re_fun(alpha_range, re_plot, tc_plot));hold on
legend('S00E35', 'S20E35', 'S40E35'); title('Static L/D');

subplot(222)
plot(alpha_range, BC_stat_S00E35_free.cl_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot)); hold on;
plot(alpha_range, BC_stat_S20E35_free.cl_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot)); hold on;
plot(alpha_range, BC_stat_S40E35_free.cl_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot)); hold on;
legend('S00E35', 'S20E35', 'S40E35'); title('Prob Cl');

subplot(224)
plot(alpha_range, BC_stat_S00E35_free.cl_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot) ./ BC_stat_S00E35_free.cd_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot) ); hold on;
plot(alpha_range, BC_stat_S20E35_free.cl_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot) ./ BC_stat_S20E35_free.cd_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot) ); hold on;
plot(alpha_range, BC_stat_S40E35_free.cl_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot) ./ BC_stat_S40E35_free.cd_of_alpha_deg_and_re_prob_fun(alpha_range, re_plot, tc_plot, mu_plot) ); hold on;
legend('S00E35', 'S20E35', 'S40E35'); title('Prob L/D');





