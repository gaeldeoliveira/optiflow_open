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
case_alfa17  = open('fig/run-over-database-alfa-17/mat/run_over_database_alfa17_intermediate_processing_20171204T045816.mat');
case_alfa20  = open('run_over_database_alfa20_gradient_20180223T062354.mat');
case_alfa21  = open('run_over_database_alfa21_gradient_20180224T024706.mat');
% Old: % ref_database = open('./data/experimental/EDB_317665fb-cd39-4321-b8cf-124f15d6f060.mat');
ref_database = open('./data/experimental/EDB_a8b44aa8-c876-4f49-8ff3-5ef3c1741218.mat');


% Define validation airfoil
% Find naca 63(4)421 in database
[i_EC_validation] = ref_database.EDB.find_cases_by_airfoil_name('NACA23021');
EC_validation_Re3 = ref_database.EDB.EC_cell{i_EC_validation(1)};
EC_validation_Re6 = ref_database.EDB.EC_cell{i_EC_validation(2)};
EC_validation_Re9 = ref_database.EDB.EC_cell{i_EC_validation(3)};
validation_airfoil_file = EC_validation_Re3.airfoil_description.airfoil_filename;


% % Receive x_theta_new, x_airfoil and parametrization orders
xgrad_alfa17 = case_alfa17.xgrad;
xgrad_alfa20 = case_alfa20.xgrad;
xgrad_alfa21 = case_alfa21.xgrad;

% % Extract airfoil polar evaluator (APE) objects for results SPR
APE_alfa21   = case_alfa21.APE;
APE_alfa20   = case_alfa20.APE;
APE_alfa17   = case_alfa17.APE;
% % Migrate corresponding system context
APE_alfa21.SC.migrate_context('a'); APE_alfa21.SW.fig_active = 0;
APE_alfa20.SC.migrate_context('b'); APE_alfa20.SW.fig_active = 0;
APE_alfa17.SC.migrate_context('c'); APE_alfa17.SW.fig_active = 0;

% % Extract airfoil polar evaluator (APE) objects for original SPR
APE_0_alfa21   = case_alfa21.GFW.APE_0;
APE_0_alfa20   = case_alfa20.GFW.APE_0;
APE_0_alfa17   = case_alfa17.GFW.APE_0;
% % Migrate corresponding system context
APE_0_alfa21.SC.migrate_context('d'); APE_0_alfa21.SW.fig_active = 0;
APE_0_alfa20.SC.migrate_context('e'); APE_0_alfa20.SW.fig_active = 0;
APE_0_alfa17.SC.migrate_context('f'); APE_0_alfa17.SW.fig_active = 0;

% % Make shape fit object for each case (they have different parametrization orders)
SF_alfa21   = shape_fit_cst(APE_alfa21.SD, 'case_alfa21');
SF_alfa20   = shape_fit_cst(APE_alfa20.SD, 'case_alfa20');
SF_alfa17   = shape_fit_cst(APE_alfa17.SD, 'case_alfa15');

% % Make shape fit object for each case (they have different parametrization orders)
x_validation_alfa21 = get_parameters_from_file(SF_alfa21 , validation_airfoil_file);
x_validation_alfa20 = get_parameters_from_file(SF_alfa20 , validation_airfoil_file);
x_validation_alfa17 = get_parameters_from_file(SF_alfa17 , validation_airfoil_file);

% Take polar conditions from some case
polar_conditions = EC_validation_Re3.polar_conditions;
%polar_conditions.N_crit = 11;

% Run polar on chosen airfoil for results SPR
[ap_alfa21   , experiment_results_alfa21  ] = APE_alfa21.make_polar_on_airfoil(  x_validation_alfa21, polar_conditions);
[ap_alfa20   , experiment_results_alfa20  ] = APE_alfa20.make_polar_on_airfoil(  x_validation_alfa20, polar_conditions);
[ap_alfa17   , experiment_results_alfa17  ] = APE_alfa17.make_polar_on_airfoil(  x_validation_alfa17, polar_conditions);

% Run polar on chosen airfoil for original SPR
[ap_0_alfa21 , experiment_results_0_alfa21] = APE_0_alfa21.make_polar_on_airfoil(x_validation_alfa21, polar_conditions);
[ap_0_alfa20 , experiment_results_0_alfa20] = APE_0_alfa21.make_polar_on_airfoil(x_validation_alfa20, polar_conditions);
[ap_0_alfa17 , experiment_results_0_alfa17] = APE_0_alfa17.make_polar_on_airfoil(x_validation_alfa17, polar_conditions);

% Show results
% Check that results are consistent for all three original closure sets
ax_ap_0 = ap_0_alfa21.plot;
          ap_0_alfa20.plot(ax_ap_0, [0.66 0.33 0.33]);
          ap_0_alfa17.plot(ax_ap_0, [0.33 0.33 0.66]);


%% Plot effect on airfoil polars

cl_al_p5 = interp1(EC_validation_Re3.al_alcl_polar_std, EC_validation_Re3.cl_alcl_polar_std, 5);
cl_al_m5 = interp1(EC_validation_Re3.al_alcl_polar_std, EC_validation_Re3.cl_alcl_polar_std,-5);
dcl_dal_deg  = (cl_al_p5 - cl_al_m5) / 10;
dcl_dal_rad  = dcl_dal_deg  * 180/pi();
dcl_dal_ratio = (2*pi()) / dcl_dal_rad;
dcl_dal_ratio2 = 1 + (dcl_dal_ratio - 1) * 0.5;

% % Plot polar at 3e6
figure(10)
hold on; grid on;
plot(EC_validation_Re3.al_alcl_polar_std / dcl_dal_ratio2 , EC_validation_Re3.cl_alcl_polar_std, 'o-k');

plot(ap_alfa21.alpha_range  , ap_alfa21.cl_alpha(ap_alfa21.alpha_range));
plot(ap_alfa20.alpha_range  , ap_alfa20.cl_alpha(ap_alfa20.alpha_range));
plot(ap_alfa17.alpha_range  , ap_alfa17.cl_alpha(ap_alfa17.alpha_range));

plot(ap_0_alfa21.alpha_range, ap_0_alfa21.cl_alpha(ap_0_alfa21.alpha_range));
%plot(ap_0_alfa20.alpha_range, ap_0_alfa20.cl_alpha(ap_0_alfa20.alpha_range));
%plot(ap_0_alfa15.alpha_range, ap_0_alfa15.cl_alpha(ap_0_alfa15.alpha_range));
%legend('Exp', 'alfa21', 'alfa20', 'alfa15', '0_alfa21', '0_alfa20', '0_alfa15');
legend('Exp', 'alfa21', 'alfa20', 'alfa17', '0_ref');

figure(11)
hold on; grid on;
plot( EC_validation_Re3.cd_clcd_polar_std          ,  EC_validation_Re3.cl_clcd_polar_std, 'o-k');
plot( ap_alfa21.cd_alpha(ap_alfa21.alpha_range)    ,  ap_alfa21.cl_alpha(ap_alfa21.alpha_range));
plot( ap_alfa20.cd_alpha(ap_alfa20.alpha_range)    ,  ap_alfa20.cl_alpha(ap_alfa20.alpha_range));
plot( ap_alfa17.cd_alpha(ap_alfa17.alpha_range)    ,  ap_alfa17.cl_alpha(ap_alfa17.alpha_range));

plot( ap_0_alfa21.cd_alpha(ap_0_alfa21.alpha_range), ap_0_alfa21.cl_alpha(ap_0_alfa21.alpha_range));
% plot( ap_0_alfa20.cd_alpha(ap_0_alfa20.alpha_range), ap_0_alfa20.cl_alpha(ap_0_alfa20.alpha_range));
% plot( ap_0_alfa15.cd_alpha(ap_0_alfa15.alpha_range), ap_0_alfa15.cl_alpha(ap_0_alfa15.alpha_range));
% legend('Exp', 'alfa21', 'alfa20', 'alfa15', '0_alfa21', '0_alfa20', '0_alfa15');
legend('Exp', 'alfa21', 'alfa20', 'alfa17', '0_ref');



%% Plot resulting closure relations

% Insert Experimental Data



rt_vinuesa = 1722  ;
hk_vinuesa = 1.74  ;
cf_vinuesa = 0.0024;

rt_vinuesa = 2255  ; %#ok<NASGU>
hk_vinuesa = 2.03  ; %#ok<NASGU>
cf_vinuesa = 0.0012; %#ok<NASGU>


% Now prepare plot of Cf closure relation
msq      = 0;
re_theta = rt_vinuesa;
hk_range = linspace(1,6);

% Compute Prelearn Reference Skin Friction Relation 
[ cf_reference] = cft(    hk_range, re_theta , msq);
% Compute Prelearn Increased Skin Friction Relation
[ cf_increased] = cft_rr( hk_range, re_theta , msq);
% Compute Postlearn Increased Skin Friction Relation
[ cf_postlearn_alfa21] = cft_ml( hk_range, re_theta , msq, case_alfa21.SPR);
[ cf_postlearn_alfa20] = cft_ml( hk_range, re_theta , msq, case_alfa20.SPR);
[ cf_postlearn_alfa17] = cft_ml( hk_range, re_theta , msq, case_alfa17.SPR);


% Plot Comparison of prelearn Cf relations
figure(87)
%subplot(221)
%plot(hk_range, cf_reference, hk_range, cf_increased, hk_range, cf_postlearn); grid on;
%plot(hk_range, cf_reference * 1e4, hk_range, cf_postlearn* 1e4, hk_vinuesa, cf_vinuesa* 1e4, 'kx'); grid on;
plot(hk_range  , cf_reference        * 1e4 , 'k-'); hold on; grid on;
plot(hk_vinuesa, cf_vinuesa          * 1e4 , 'kx'); 
plot(hk_range  , cf_postlearn_alfa21 * 1e4);
plot(hk_range  , cf_postlearn_alfa20 * 1e4);
plot(hk_range  , cf_postlearn_alfa17 * 1e4);
legend('Previous Closure', 'DNS [Vinuesa et al.]', 'Learned Closure alfa21', 'Learned Closure alfa20', 'Learned Closure alfa15');

xlabel('H_k - Kinematic Shape Factor');
ylabel('C_f x 10^4');
title(['Skin Friction (Re_{\theta} = ' num2str(re_theta) , ')']);
axis([1 5 -20 80])

% % Now go on to explore H* correlation

% Compute Prelearn Reference H* Skin Friction Relation 
[ hs_reference] = hst(    hk_range, re_theta , msq);
% Compute Prelearn Increased H* Skin Friction Relation
[ hs_increased] = hstivw( hk_range, re_theta , msq);
% Compute Postlearn Increased H* Skin Friction Relation
%[ hs_postlearn] = hst_ml( hk_range, re_theta , msq, CS1.SPR);
[ hs_postlearn_alfa21] = hst_ml( hk_range, re_theta , msq, case_alfa21.SPR);
[ hs_postlearn_alfa20] = hst_ml( hk_range, re_theta , msq, case_alfa20.SPR);
[ hs_postlearn_alfa17] = hst_ml( hk_range, re_theta , msq, case_alfa17.SPR);

% Plot Comparison of prelearn Cf relations
 figure(85)
%subplot(222)
%plot(hk_range, hs_reference, hk_range, hs_increased, hk_range, hs_postlearn); grid on;
%plot(hk_range, hs_reference, hk_range, hs_postlearn); grid on;
plot(hk_range  , hs_reference        , 'k-'); hold on; grid on;
plot(hk_range  , hs_postlearn_alfa21 );
plot(hk_range  , hs_postlearn_alfa20 );
plot(hk_range  , hs_postlearn_alfa17 );
legend('Previous Closure', 'Learned Closure alfa21', 'Learned Closure alfa20', 'Learned Closure alfa15');
xlabel('H_k - Kinematic Shape Factor');
ylabel('H^{*} - Skin Friction Coefficient');
title(['Energy Shape Factor (Re_{\theta} = ' , num2str(re_theta), ')']);
axis([1 5 1.0 2])

du96 = [100,0;99.5215000000000,0.120700000000000;98.9140000000000,0.261600000000000;98.1488000000000,0.439600000000000;97.2238000000000,0.655500000000000;96.1816000000000,0.898600000000000;95.0733000000000,1.15560000000000;93.9321000000000,1.41900000000000;92.7768000000000,1.68540000000000;91.6166000000000,1.95260000000000;90.4552000000000,2.21940000000000;89.2932000000000,2.48540000000000;88.1321000000000,2.75040000000000;86.9710000000000,3.01380000000000;85.8102000000000,3.27610000000000;84.6502000000000,3.53630000000000;83.4886000000000,3.79510000000000;82.3266000000000,4.05290000000000;81.1653000000000,4.30930000000000;80.0050000000000,4.56360000000000;78.8446000000000,4.81580000000000;77.6840000000000,5.06610000000000;76.5221000000000,5.31430000000000;75.3593000000000,5.56140000000000;74.1966000000000,5.80660000000000;73.0338000000000,6.05000000000000;71.8699000000000,6.29160000000000;70.7048000000000,6.53200000000000;69.5396000000000,6.77090000000000;68.3737000000000,7.00850000000000;67.2084000000000,7.24440000000000;66.0420000000000,7.47880000000000;64.8756000000000,7.71240000000000;63.7112000000000,7.94390000000000;62.5483000000000,8.17300000000000;61.3873000000000,8.39890000000000;60.2270000000000,8.62140000000000;59.0675000000000,8.84050000000000;57.9099000000000,9.05560000000000;56.7542000000000,9.26600000000000;55.6008000000000,9.47110000000000;54.4494000000000,9.67020000000000;53.3005000000000,9.86270000000000;52.1544000000000,10.0478000000000;51.0115000000000,10.2244000000000;49.8708000000000,10.3919000000000;48.7327000000000,10.5497000000000;47.5977000000000,10.6966000000000;46.4634000000000,10.8320000000000;45.3285000000000,10.9564000000000;44.1961000000000,11.0696000000000;43.0660000000000,11.1697000000000;41.9368000000000,11.2572000000000;40.8093000000000,11.3312000000000;39.6833000000000,11.3916000000000;38.5591000000000,11.4375000000000;37.4357000000000,11.4691000000000;36.3150000000000,11.4859000000000;35.1964000000000,11.4871000000000;34.0792000000000,11.4726000000000;32.9641000000000,11.4421000000000;31.8499000000000,11.3952000000000;30.7361000000000,11.3320000000000;29.6225000000000,11.2528000000000;28.5106000000000,11.1579000000000;27.3999000000000,11.0467000000000;26.2903000000000,10.9199000000000;25.1826000000000,10.7771000000000;24.0762000000000,10.6187000000000;22.9731000000000,10.4450000000000;21.8727000000000,10.2556000000000;20.7762000000000,10.0509000000000;19.6846000000000,9.83040000000000;18.5974000000000,9.59410000000000;17.5158000000000,9.34190000000000;16.4406000000000,9.07370000000000;15.3732000000000,8.78940000000000;14.3144000000000,8.48830000000000;13.2650000000000,8.17030000000000;12.2261000000000,7.83470000000000;11.1990000000000,7.48150000000000;10.1860000000000,7.11030000000000;9.18920000000000,6.72100000000000;8.21140000000000,6.31300000000000;7.25620000000000,5.88680000000000;6.32880000000000,5.44280000000000;5.43520000000000,4.98220000000000;4.58510000000000,4.50870000000000;3.79090000000000,4.02810000000000;3.06700000000000,3.54950000000000;2.42820000000000,3.08600000000000;1.88360000000000,2.64960000000000;1.43300000000000,2.24870000000000;1.06880000000000,1.88700000000000;0.778300000000000,1.56260000000000;0.548200000000000,1.27060000000000;0.367300000000000,1.00600000000000;0.227500000000000,0.764100000000000;0.122800000000000,0.540500000000000;0.0501000000000000,0.332400000000000;0.00900000000000000,0.136200000000000;0.00110000000000000,-0.0466000000000000;0.0306000000000000,-0.232900000000000;0.106200000000000,-0.421400000000000;0.224200000000000,-0.611000000000000;0.381000000000000,-0.804700000000000;0.578900000000000,-1.00610000000000;0.823500000000000,-1.21870000000000;1.12440000000000,-1.44600000000000;1.49420000000000,-1.69140000000000;1.94790000000000,-1.95720000000000;2.49860000000000,-2.24250000000000;3.15440000000000,-2.54280000000000;3.91230000000000,-2.85100000000000;4.75800000000000,-3.15700000000000;5.67360000000000,-3.45260000000000;6.64170000000000,-3.73310000000000;7.64810000000000,-3.99610000000000;8.68280000000000,-4.24090000000000;9.73940000000000,-4.46770000000000;10.8148000000000,-4.67780000000000;11.9036000000000,-4.87250000000000;13.0033000000000,-5.05190000000000;14.1130000000000,-5.21750000000000;15.2317000000000,-5.37000000000000;16.3574000000000,-5.51080000000000;17.4888000000000,-5.64010000000000;18.6252000000000,-5.75870000000000;19.7662000000000,-5.86710000000000;20.9113000000000,-5.96590000000000;22.0603000000000,-6.05580000000000;23.2114000000000,-6.13730000000000;24.3638000000000,-6.21000000000000;25.5184000000000,-6.27440000000000;26.6747000000000,-6.33080000000000;27.8319000000000,-6.37940000000000;28.9893000000000,-6.42020000000000;30.1475000000000,-6.45290000000000;31.3069000000000,-6.47790000000000;32.4666000000000,-6.49540000000000;33.6264000000000,-6.50540000000000;34.7859000000000,-6.50780000000000;35.9451000000000,-6.50240000000000;37.1048000000000,-6.48920000000000;38.2646000000000,-6.46860000000000;39.4241000000000,-6.44020000000000;40.5840000000000,-6.40420000000000;41.7452000000000,-6.36070000000000;42.9071000000000,-6.31040000000000;44.0687000000000,-6.25310000000000;45.2303000000000,-6.18880000000000;46.3922000000000,-6.11770000000000;47.5537000000000,-6.04000000000000;48.7147000000000,-5.95530000000000;49.8769000000000,-5.86380000000000;51.0389000000000,-5.76630000000000;52.2003000000000,-5.66230000000000;53.3625000000000,-5.55190000000000;54.5242000000000,-5.43570000000000;55.6855000000000,-5.31320000000000;56.8470000000000,-5.18480000000000;58.0092000000000,-5.05040000000000;59.1721000000000,-4.91060000000000;60.3352000000000,-4.76540000000000;61.4989000000000,-4.61520000000000;62.6626000000000,-4.46020000000000;63.8268000000000,-4.30030000000000;64.9924000000000,-4.13580000000000;66.1589000000000,-3.96720000000000;67.3273000000000,-3.79460000000000;68.4974000000000,-3.61870000000000;69.6700000000000,-3.43970000000000;70.8454000000000,-3.25860000000000;72.0224000000000,-3.07600000000000;73.2001000000000,-2.89270000000000;74.3763000000000,-2.70910000000000;75.5515000000000,-2.52690000000000;76.7255000000000,-2.34640000000000;77.8970000000000,-2.16830000000000;79.0664000000000,-1.99350000000000;80.2352000000000,-1.82200000000000;81.4021000000000,-1.65400000000000;82.5663000000000,-1.49060000000000;83.7292000000000,-1.33200000000000;84.8898000000000,-1.17870000000000;86.0490000000000,-1.03160000000000;87.2068000000000,-0.890500000000000;88.3621000000000,-0.756300000000000;89.5137000000000,-0.629600000000000;90.6625000000000,-0.511600000000000;91.8067000000000,-0.402500000000000;92.9432000000000,-0.304000000000000;94.0697000000000,-0.217200000000000;95.1763000000000,-0.143400000000000;96.2449000000000,-0.0848000000000000;97.2468000000000,-0.0435000000000000;98.1424000000000,-0.0185000000000000;98.8961000000000,-0.00840000000000000;99.5085000000000,-0.00690000000000000;100,0]
save(



