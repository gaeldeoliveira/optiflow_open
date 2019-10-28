%x_results_alfa21 =  [0.9987    0.9436    0.9105    1.1407    0.9984    0.9584    1.0207    0.9839]


% Add Necessary Paths
addpath src
addpath user_src
addpath dev
addpath dev/0_closure_relations/

% Set name of validation airfoil
% validation_airfoil_file = 'data/experimental/verification/du91.txt';
% validation_airfoil_file = './data/experimental/NACA/TR824-Airfoils/naca23021.air';
% validation_airfoil_file = 'data/experimental/verification/naca634421_coord.txt';

% Load results from a gradient solution
%case_alfa15  = open('fig/run-over-database-alfa-15/mat/run_over_database_alfa15_20171123T105445.mat');
%case_alfa17  = open('fig/run-over-database-alfa-17/mat/run_over_database_alfa17_intermediate_processing_20171204T045816.mat');
%case_alfa20  = open('run_over_database_alfa20_gradient_20180223T062354.mat');
case_alfa21  = open('run_over_database_alfa21_gradient_20180224T024706.mat');
% Old: % ref_database = open('./data/experimental/EDB_317665fb-cd39-4321-b8cf-124f15d6f060.mat');
ref_database = open('./data/experimental/EDB_a8b44aa8-c876-4f49-8ff3-5ef3c1741218.mat');

% Define validation airfoil
[i_NACA2415] = ref_database.EDB.find_cases_by_airfoil_name('NACA2415');
EC_NACA2415_9e6 = ref_database.EDB.EC_cell{i_NACA2415(3)};

% Define validation airfoil
validation_airfoil_file = EC_NACA2415_9e6.airfoil_description.airfoil_filename;

% % Receive x_theta_new, x_airfoil and parametrization orders
xgrad_alfa21 = case_alfa21.xgrad;

% % Extract airfoil polar evaluator (APE) objects for results SPR
APE_alfa21   = case_alfa21.APE;
% % Migrate corresponding system context
APE_alfa21.SC.N_cores   = 1; APE_alfa21.SC.migrate_context('a'); APE_alfa21.SW.fig_active = 0;

% % Extract airfoil polar evaluator (APE) objects for original SPR
APE_0_alfa21   = case_alfa21.GFW.APE_0;
% % Migrate corresponding system context
APE_0_alfa21.SC.N_cores = 1; APE_0_alfa21.SC.migrate_context('d'); APE_0_alfa21.SW.fig_active = 0;

% % Make shape fit object for each case (they have different parametrization orders)
SF_alfa21   = shape_fit_cst(APE_alfa21.SD, 'case_alfa21');

% % Make shape fit object for each case (they have different parametrization orders)
x_validation_alfa21 = get_parameters_from_file(SF_alfa21 , validation_airfoil_file);

% Define conditions of experimental reference
EPC                 = EC_NACA2415_9e6.polar_conditions;
EPC.INGE            = false;    % Van ingen method :)
polar_conditions    = EPC;
EPC.INGE            = true;    % Van ingen method :)
polar_conditions_eN = EPC;

% Run polar on chosen airfoil for results SPR
[ap_alfa21     , experiment_results_alfa21      ] = APE_alfa21.make_polar_on_airfoil(  x_validation_alfa21, polar_conditions   );
[ap_alfa21_eN  , experiment_results_alfa21_eN   ] = APE_alfa21.make_polar_on_airfoil(  x_validation_alfa21, polar_conditions_eN);
%[ap_alfa20   , experiment_results_alfa20  ] = APE_alfa20.make_polar_on_airfoil(  x_validation_alfa20, polar_conditions);
%[ap_alfa17   , experiment_results_alfa17  ] = APE_alfa17.make_polar_on_airfoil(  x_validation_alfa17, polar_conditions);

% Run polar on chosen airfoil for original SPR
[ap_0_alfa21    , experiment_results_0_alfa21   ] = APE_0_alfa21.make_polar_on_airfoil(x_validation_alfa21, polar_conditions   );
[ap_0_alfa21_eN , experiment_results_0_alfa21_eN] = APE_0_alfa21.make_polar_on_airfoil(x_validation_alfa21, polar_conditions_eN);
%[ap_0_alfa20 , experiment_results_0_alfa20] = APE_0_alfa20.make_polar_on_airfoil(x_validation_alfa20, polar_conditions);
%[ap_0_alfa17 , experiment_results_0_alfa17] = APE_0_alfa17.make_polar_on_airfoil(x_validation_alfa17, polar_conditions);

% Show results
% Check that results are consistent for all three original closure sets
ax_ap_0 = ap_0_alfa21.plot;
          ap_0_alfa21_eN.plot(ax_ap_0, [0.66 0.33 0.33])
            ap_alfa21.plot(ax_ap_0, [0.33 0.33 0.66]);
            ap_alfa21_eN.plot(ax_ap_0, [0.33 0.33 0.66]);
            

% Conclusion: will kick out alfa20;

%% Plot effect on airfoil polars

% % Plot polar at 3e6
figure(10)
hold on; grid on;
plot(EC_NACA2415_9e6.al_alcl_polar_std, EC_NACA2415_9e6.cl_alcl_polar_std, 'o-k');

plot(ap_alfa21.alpha_range  , ap_alfa21.cl_alpha(ap_alfa21.alpha_range));
plot(ap_alfa21_eN.alpha_range  , ap_alfa21_eN.cl_alpha(ap_alfa21_eN.alpha_range));
%plot(ap_alfa20.alpha_range  , ap_alfa20.cl_alpha(ap_alfa20.alpha_range));
%plot(ap_alfa17.alpha_range  , ap_alfa17.cl_alpha(ap_alfa17.alpha_range));

plot(ap_0_alfa21.alpha_range, ap_0_alfa21.cl_alpha(ap_0_alfa21.alpha_range));
plot(ap_0_alfa21_eN.alpha_range, ap_0_alfa21_eN.cl_alpha(ap_0_alfa21_eN.alpha_range));
%plot(ap_0_alfa20.alpha_range, ap_0_alfa20.cl_alpha(ap_0_alfa20.alpha_range));
%plot(ap_0_alfa15.alpha_range, ap_0_alfa15.cl_alpha(ap_0_alfa15.alpha_range));
%legend('Exp', 'alfa21', 'alfa20', 'alfa15', '0_alfa21', '0_alfa20', '0_alfa15');
%legend('Exp', 'alfa21', 'alfa20', 'alfa17', '0_ref');
%legend('Exp', 'alfa21', '0_ref');
legend('Exp', 'alfa21', 'alfa21_eN', '0_ref', '0_ref_eN');

figure(11)
hold on; grid on;
plot(EC_NACA2415_9e6.cd_alclcd_polar_std, EC_NACA2415_9e6.cl_alclcd_polar_std, 'o-k');

plot( ap_alfa21.cd_alpha(ap_alfa21.alpha_range)    ,  ap_alfa21.cl_alpha(ap_alfa21.alpha_range));
plot( ap_alfa21_eN.cd_alpha(ap_alfa21_eN.alpha_range)    ,  ap_alfa21_eN.cl_alpha(ap_alfa21_eN.alpha_range));

plot( ap_0_alfa21.cd_alpha(ap_0_alfa21.alpha_range), ap_0_alfa21.cl_alpha(ap_0_alfa21.alpha_range));
plot( ap_0_alfa21_eN.cd_alpha(ap_0_alfa21_eN.alpha_range), ap_0_alfa21_eN.cl_alpha(ap_0_alfa21_eN.alpha_range));
%legend('Exp', 'alfa21', 'alfa20', 'alfa17', '0_ref');
%legend('Exp', 'alfa21', '0_ref');
legend('Exp', 'alfa21', 'alfa21_eN', '0_ref', '0_ref_eN');

% Conclusion: will prefer alfa21 to alfa17
% 2nd Conclusion: chosen nonlinear eN method
%
%
%   Mention this: non-linear eN method poses convergence challenges for
%   gradient minimization, and the learning process was conducted . Polar
%   curves display eh
%
%

% % Now prepare final polar plots (use non-linear eN where available, fallback to linear eN other wise , for lift only)
% Make composite for cl of ap_alfa21 case
ap_alfa21_alpha_lin   = ap_alfa21.alpha_range;
ap_alfa21_alpha_eN    = ap_alfa21_eN.alpha_range;
ap_alfa21_cl_lin      = ap_alfa21.cl_alpha(   ap_alfa21.alpha_range);
ap_alfa21_cl_eN       = ap_alfa21_eN.cl_alpha(ap_alfa21_alpha_eN   );

i_below               = ap_alfa21_alpha_lin < min(ap_alfa21_alpha_eN);
i_above               = ap_alfa21_alpha_lin > max(ap_alfa21_alpha_eN);

ap_alfa21_alpha_mixed = [ap_alfa21_alpha_lin(i_below) , ap_alfa21_alpha_eN, ap_alfa21_alpha_lin(i_above)];
ap_alfa21_cl_mixed    = [ap_alfa21_cl_lin(i_below) , ap_alfa21_cl_eN, ap_alfa21_cl_lin(i_above)];

% Make composite for cl of ap_0_alfa21 case
ap_0_alfa21_alpha_lin   = ap_0_alfa21.alpha_range;
ap_0_alfa21_alpha_eN    = ap_0_alfa21_eN.alpha_range;
ap_0_alfa21_cl_lin      = ap_0_alfa21.cl_alpha(   ap_0_alfa21.alpha_range);
ap_0_alfa21_cl_eN       = ap_0_alfa21_eN.cl_alpha(ap_0_alfa21_alpha_eN   );

i_below               = ap_0_alfa21_alpha_lin < min(ap_0_alfa21_alpha_eN);
i_above               = ap_0_alfa21_alpha_lin > max(ap_0_alfa21_alpha_eN);

ap_0_alfa21_alpha_mixed = [ap_0_alfa21_alpha_lin(i_below) , ap_0_alfa21_alpha_eN, ap_0_alfa21_alpha_lin(i_above)];
ap_0_alfa21_cl_mixed    = [ap_0_alfa21_cl_lin(i_below) , ap_0_alfa21_cl_eN, ap_0_alfa21_cl_lin(i_above)];

% % Now test with plots
figure(110)
hold on; grid on;
plot(EC_NACA2415_9e6.al_alcl_polar_std, EC_NACA2415_9e6.cl_alcl_polar_std, 'o-k');

plot(ap_alfa21.alpha_range  , ap_alfa21.cl_alpha(ap_alfa21.alpha_range));
plot(ap_alfa21_eN.alpha_range  , ap_alfa21_eN.cl_alpha(ap_alfa21_eN.alpha_range));
plot(ap_alfa21_alpha_mixed , ap_alfa21_cl_mixed, '--');
% Looks OK!

figure(210)
hold on; grid on;
plot(EC_NACA2415_9e6.al_alcl_polar_std, EC_NACA2415_9e6.cl_alcl_polar_std, 'o-k');

plot(ap_0_alfa21.alpha_range  , ap_0_alfa21.cl_alpha(ap_0_alfa21.alpha_range));
plot(ap_0_alfa21_eN.alpha_range  , ap_0_alfa21_eN.cl_alpha(ap_0_alfa21_eN.alpha_range));
plot(ap_0_alfa21_alpha_mixed , ap_0_alfa21_cl_mixed, '--');
% Looks OK!

% % Now draft of final plot of DU96W180 polar
figure(310)
hold on; grid on;
plot(EC_NACA2415_9e6.al_alcl_polar_std, EC_NACA2415_9e6.cl_alcl_polar_std, 'o-k');
plot(ap_0_alfa21_alpha_mixed , ap_0_alfa21_cl_mixed, '--');
plot(ap_alfa21_alpha_mixed , ap_alfa21_cl_mixed, '--');
legend('Experimental', 'Original RFOIL', 'Learned Rfoil');

figure(311)
hold on; grid on;
plot(EC_NACA2415_9e6.cd_alclcd_polar_std, EC_NACA2415_9e6.cl_alclcd_polar_std, 'o-k');
plot( ap_0_alfa21_eN.cd_alpha(ap_0_alfa21_eN.alpha_range), ap_0_alfa21_eN.cl_alpha(ap_0_alfa21_eN.alpha_range));
plot( ap_alfa21_eN.cd_alpha(ap_alfa21_eN.alpha_range)    ,  ap_alfa21_eN.cl_alpha(ap_alfa21_eN.alpha_range));
legend('Experimental', 'Original RFOIL', 'Learned Rfoil');

% % Now make it into a side-by-side polar curve
figure(312)
%subplot(122)
subplot(4,2,[2 4 6])
hold on; grid on;
plot(EC_NACA2415_9e6.al_alcl_polar_std, EC_NACA2415_9e6.cl_alcl_polar_std, 'o-k');
plot(ap_0_alfa21_alpha_mixed , ap_0_alfa21_cl_mixed);
plot(ap_alfa21_alpha_mixed , ap_alfa21_cl_mixed);
legend('Langley LTPT', 'Original Rfoil', 'Learned Rfoil', 'Location', 'NorthWest');
axis([-5 20 -0.4 2.0]);
title('Lift Polar - NACA2415 @ Re =9e6');
xlabel('\alpha (deg.)');
ylabel('C_l');
%subplot(121)
subplot(4,2,[1 3 5])
hold on; grid on;
plot(EC_NACA2415_9e6.cd_alclcd_polar_std * 10^4, EC_NACA2415_9e6.cl_alclcd_polar_std, 'o-k');
plot( ap_0_alfa21_eN.cd_alpha(ap_0_alfa21_eN.alpha_range) * 10^4, ap_0_alfa21_eN.cl_alpha(ap_0_alfa21_eN.alpha_range));
plot( ap_alfa21_eN.cd_alpha(ap_alfa21_eN.alpha_range)     * 10^4,  ap_alfa21_eN.cl_alpha(ap_alfa21_eN.alpha_range));
legend('Langley LTPT', 'Original Rfoil', 'Learned Rfoil', 'Location', 'NorthWest');
axis([0.0020*10^4 0.0180*10^4 -0.4 2.0]);% axis([0.00*10^4 0.015*10^4 -0.4 2.0]);
title('Drag Polar - NACA2415 @ Re =9e6');
xlabel('C_d x 10^4');
ylabel('C_l');

% % And write figure to file
set(gcf, 'PaperType', 'A5');
orient landscape;
print -dpdf -r300 torque2018ML_NACA2415_polar.pdf






