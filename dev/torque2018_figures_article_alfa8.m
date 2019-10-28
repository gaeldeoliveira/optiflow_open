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

% % Define conditions of experimental reference
% EPC                 = EC_NACA2415_9e6.polar_conditions;
% EPC.INGE            = false;    % Van ingen method :)
% polar_conditions    = EPC;
% EPC.INGE            = true;    % Van ingen method :)
% polar_conditions_eN = EPC;
% 
% % Run polar on chosen airfoil for results SPR
% [ap_alfa21     , experiment_results_alfa21      ] = APE_alfa21.make_polar_on_airfoil(  x_validation_alfa21, polar_conditions   );
% [ap_alfa21_eN  , experiment_results_alfa21_eN   ] = APE_alfa21.make_polar_on_airfoil(  x_validation_alfa21, polar_conditions_eN);
% %[ap_alfa20   , experiment_results_alfa20  ] = APE_alfa20.make_polar_on_airfoil(  x_validation_alfa20, polar_conditions);
% %[ap_alfa17   , experiment_results_alfa17  ] = APE_alfa17.make_polar_on_airfoil(  x_validation_alfa17, polar_conditions);
% 
% % Run polar on chosen airfoil for original SPR
% [ap_0_alfa21    , experiment_results_0_alfa21   ] = APE_0_alfa21.make_polar_on_airfoil(x_validation_alfa21, polar_conditions   );
% [ap_0_alfa21_eN , experiment_results_0_alfa21_eN] = APE_0_alfa21.make_polar_on_airfoil(x_validation_alfa21, polar_conditions_eN);
% %[ap_0_alfa20 , experiment_results_0_alfa20] = APE_0_alfa20.make_polar_on_airfoil(x_validation_alfa20, polar_conditions);
% %[ap_0_alfa17 , experiment_results_0_alfa17] = APE_0_alfa17.make_polar_on_airfoil(x_validation_alfa17, polar_conditions);
% 
% % Show results
% % Check that results are consistent for all three original closure sets
% ax_ap_0 = ap_0_alfa21.plot;
%           ap_0_alfa21_eN.plot(ax_ap_0, [0.66 0.33 0.33])
%             ap_alfa21.plot(ax_ap_0, [0.33 0.33 0.66]);
%             ap_alfa21_eN.plot(ax_ap_0, [0.33 0.33 0.66]);
%             
%
figure(54)
subplot(221)
ClBinEdges        = [0.0200 0.0400 0.0600 0.0800 0.1000 0.1200];
H_APE_0_alfa21_cl = histogram(APE_0_alfa21.cl_accuracy_metric_array, ClBinEdges); hold on;
H_APE_alfa21_cl   = histogram(APE_alfa21.cl_accuracy_metric_array  , ClBinEdges);
title(['Global Cl accuracy =' num2str(case_alfa21.normalized_accuracy_profile.cl_global_accuracy, 3)]);
xlabel('RMS error of lift coefficient'); ylabel('# of Polar Curves');
H_APE_alfa21_cl.Parent.YLim = [0 14];
grid on;
legend('Before learning' , 'After learning')
% legend('Before learning' , 'After learning', 'Location', 'South')

subplot(222)
% CdBinEdges        = [0 3.0000e-04 6.0000e-04 9.0000e-04 0.0012 0.0015 0.0018];
% %CdBinEdges        = [0 5.0000e-04 1.0000e-03 0.0015]
% CdBinEdges        = [0 3.5000e-04 7.0000e-04 0.0010 0.0014 0.0017];
CdBinEdges        = [2.0000e-04 4.6000e-04 7.2000e-04 9.8000e-04 0.0012 0.0015] * 10^4
H_APE_0_alfa21_cd = histogram(APE_0_alfa21.cd_accuracy_metric_array * 10^4, CdBinEdges); hold on;
H_APE_alfa21_cd   = histogram(APE_alfa21.cd_accuracy_metric_array   * 10^4, CdBinEdges);
title(['Global Cd accuracy =' num2str(case_alfa21.normalized_accuracy_profile.cd_global_accuracy, 3)]);
xlabel('RMS error of drag coefficient x 10^4'); ylabel('# of polar curves');
H_APE_alfa21_cd.Parent.YLim = [0 14];
H_APE_alfa21_cd.Parent.XTick = CdBinEdges;
grid on;
legend('Before learning' , 'After learning');
% legend('Before learning' , 'After learning', 'Location', 'South')

% % Write to figure!
set(gcf, 'PaperType', 'A5');
orient landscape;
print -dpdf -r300 torque2018ML_histograms.pdf

% Edit legend positions
% % % Write to figure!
% set(gcf, 'PaperType', 'A5');
% orient landscape;
% print -dpdf -r300 torque2018ML_histograms3.pdf



% Conclusion: will kick out alfa20;

%% Plot effect on airfoil polars

% % 






