addpath src

ICM = struct()                   ;      % make struct
ICM.C_WT_fix     =      14272366.8 + 0.8 * 16417084.0  + 14272366.8;      % [Eur]
ICM.C_WT_var     =  0.2*16417084.0 + 0.0 * 14272366.8  ;      % [Eur] Variable cost component
ICM.C_WT_fix     =      14272366.8 + 0.4 * 16417084.0  + 14272366.8;      % [Eur]
ICM.C_WT_var     =  0.4*16417084.0 + 0.0 * 14272366.8  ;      % [Eur] Variable cost component
ICM.C_WT_fix     =      14272366.8 + 0.7 * 16417084.0  ;      % [Eur]
ICM.C_WT_var     =  0.3*16417084.0 + 0.0 * 14272366.8  ;      % [Eur] Variable cost component

ICM.LCOM         = 31.09;                                        % [Eur/MWh] Levelized cost of maintenance


ICM.gamma_RBM    =        1.00   ;      % Exponent of root bending moment scaling

ICM.r_yr         =         0.0738;      % Yearly interest rate of cost of capital (7.38% nominal, 2% inflation, 5.38% real in Innwind model) (nominal leads to right LCOE (62.09 is spreadsheet, 61.6 here) (pre-maintenace)
ICM.N_yr         =        25     ;      % Number of years turbine operates

ICM.a_low_cp     =        25068 ;       % [MWh] 
ICM.b_low_cp     =        40482 ;       % [MWh] 

ICM.a_high_cp    =        31196 ;       % [MWh] 
ICM.b_high_cp    =        28774 ;       % [MWh] 

ICM.CP_threshold =         0.5234;      % [adim.]

ICM.m0           =         41716;       % [Kg]

% Load reference case to find original design root bending moment
date_str         = 'DO_20180915T200345';
ref_case_ICM     = load([date_str '/prob_FFA_free.mat'   ]);
cf_ref_design    = ref_case_ICM.cf_wrapper_fun(ref_case_ICM.A_free_0);
ICM.CQ_RBM_0     =  cf_ref_design(2);
ICM.CP_0         = -cf_ref_design(1);


% Make discount factor
n_year_vector        = 1:ICM.N_yr;
discount_year_vector = (1+ICM.r_yr).^(-n_year_vector);
ICM.discounted_sum   = sum(discount_year_vector) ;%/ ICM.N_yr;

% Expose AEP function to outer world
AEP_fun              = @(CP)         AEP_internal_fun(ICM, CP);

% Make function for total cost of wind turbine 
C_WT_total_fun       = @(CQ_RBM)     ICM.C_WT_fix + ICM.C_WT_var * (CQ_RBM ./ ICM.CQ_RBM_0).^ICM.gamma_RBM;
% Make function for levelized cost of energy (no maintenance costs)
LCOE_fun             = @(CQ_RBM, CP) C_WT_total_fun(CQ_RBM) / (AEP_fun(CP) * ICM.discounted_sum) + ICM.LCOM;
% Store cost of energy of original design
ICM.LCOE_0           = LCOE_fun(ICM.CQ_RBM_0, ICM.CP_0);
disp('LCOE_0:')
disp(ICM.LCOE_0)

% Now make function for CQ_RBM given LCOE and CP_rotor
CQ_RBM_iso_cost_fun  = @(LCOE, CP) ICM.CQ_RBM_0 * ( ((LCOE-ICM.LCOM)./ICM.C_WT_var) .* AEP_fun(CP) .* ICM.discounted_sum -  ICM.C_WT_fix/ICM.C_WT_var ).^(1/ICM.gamma_RBM);

% Now verify consistency of iso-cost line with LCOE function 
% (PASS corresponds to true (1), given absence of rounding error! check
% manually otherwise!)
disp('IsoCost inconsistency residual:')
disp((CQ_RBM_iso_cost_fun(ICM.LCOE_0, ICM.CP_0) - ICM.CQ_RBM_0) / ICM.CQ_RBM_0 * 100)

% % And now plot!
% Make a range of Cp values
CP_isocost_range = linspace(0.3, 16/27);

LCOE_m10 = ICM.LCOE_0 * (1 - 0.040);
LCOE_m09 = ICM.LCOE_0 * (1 - 0.036);
LCOE_m08 = ICM.LCOE_0 * (1 - 0.032);
LCOE_m07 = ICM.LCOE_0 * (1 - 0.028);
LCOE_m06 = ICM.LCOE_0 * (1 - 0.024);
LCOE_m05 = ICM.LCOE_0 * (1 - 0.020);
LCOE_m04 = ICM.LCOE_0 * (1 - 0.016);
LCOE_m03 = ICM.LCOE_0 * (1 - 0.012);
LCOE_m02 = ICM.LCOE_0 * (1 - 0.008);
LCOE_m01 = ICM.LCOE_0 * (1 - 0.004);
LCOE_000 = ICM.LCOE_0 ;
LCOE_p01 = ICM.LCOE_0 * (1 + 0.004);
LCOE_p02 = ICM.LCOE_0 * (1 + 0.008);
LCOE_p03 = ICM.LCOE_0 * (1 + 0.012);
LCOE_p04 = ICM.LCOE_0 * (1 + 0.016);
LCOE_p05 = ICM.LCOE_0 * (1 + 0.020);
LCOE_p06 = ICM.LCOE_0 * (1 + 0.024);
LCOE_p07 = ICM.LCOE_0 * (1 + 0.028);
LCOE_p08 = ICM.LCOE_0 * (1 + 0.032);
LCOE_p09 = ICM.LCOE_0 * (1 + 0.036);
LCOE_p10 = ICM.LCOE_0 * (1 + 0.040);


CQ_RBM_iso_cost_m10 = CQ_RBM_iso_cost_fun(LCOE_m10, CP_isocost_range);
CQ_RBM_iso_cost_m09 = CQ_RBM_iso_cost_fun(LCOE_m09, CP_isocost_range);
CQ_RBM_iso_cost_m08 = CQ_RBM_iso_cost_fun(LCOE_m08, CP_isocost_range);
CQ_RBM_iso_cost_m07 = CQ_RBM_iso_cost_fun(LCOE_m07, CP_isocost_range);
CQ_RBM_iso_cost_m06 = CQ_RBM_iso_cost_fun(LCOE_m06, CP_isocost_range);
CQ_RBM_iso_cost_m05 = CQ_RBM_iso_cost_fun(LCOE_m05, CP_isocost_range);
CQ_RBM_iso_cost_m04 = CQ_RBM_iso_cost_fun(LCOE_m04, CP_isocost_range);
CQ_RBM_iso_cost_m03 = CQ_RBM_iso_cost_fun(LCOE_m03, CP_isocost_range);
CQ_RBM_iso_cost_m02 = CQ_RBM_iso_cost_fun(LCOE_m02, CP_isocost_range);
CQ_RBM_iso_cost_m01 = CQ_RBM_iso_cost_fun(LCOE_m01, CP_isocost_range);
CQ_RBM_iso_cost_000 = CQ_RBM_iso_cost_fun(LCOE_000, CP_isocost_range);
CQ_RBM_iso_cost_p01 = CQ_RBM_iso_cost_fun(LCOE_p01, CP_isocost_range);
CQ_RBM_iso_cost_p02 = CQ_RBM_iso_cost_fun(LCOE_p02, CP_isocost_range);
CQ_RBM_iso_cost_p03 = CQ_RBM_iso_cost_fun(LCOE_p03, CP_isocost_range);
CQ_RBM_iso_cost_p04 = CQ_RBM_iso_cost_fun(LCOE_p04, CP_isocost_range);
CQ_RBM_iso_cost_p05 = CQ_RBM_iso_cost_fun(LCOE_p05, CP_isocost_range);
CQ_RBM_iso_cost_p06 = CQ_RBM_iso_cost_fun(LCOE_p06, CP_isocost_range);
CQ_RBM_iso_cost_p07 = CQ_RBM_iso_cost_fun(LCOE_p07, CP_isocost_range);
CQ_RBM_iso_cost_p08 = CQ_RBM_iso_cost_fun(LCOE_p08, CP_isocost_range);
CQ_RBM_iso_cost_p09 = CQ_RBM_iso_cost_fun(LCOE_p09, CP_isocost_range);
CQ_RBM_iso_cost_p10 = CQ_RBM_iso_cost_fun(LCOE_p10, CP_isocost_range);

% LCOE_isocost_stages = ICM.LCOE_0 * (1 + [-0.032:0.004:0.032]);
% CQ_RBM_iso_cost     = 


% Prepare to overlay plots
% Get data from 

date_str = 'DO_20180915T200345';

prob_FFA_free            = load([date_str '/prob_FFA_free.mat'           ]);
prob_high_glide_free     = load([date_str '/prob_high_glide_free.mat'    ]);
prob_smooth_glide_free   = load([date_str '/prob_smooth_glide_free.mat'  ]);
prob_free_limit_case     = load([date_str '/prob_FFA_free_limit_case.mat']);


% Reorder pareto fronts (prob)
[prob_FFA_free.fval_CP_sorted          , prob_FFA_free.sort_index         ] = sort(prob_FFA_free.fval(         :,1)); prob_FFA_free.fval_CQ_sorted          = prob_FFA_free.fval(         prob_FFA_free.sort_index         ,2);
[prob_high_glide_free.fval_CP_sorted   , prob_high_glide_free.sort_index  ] = sort(prob_high_glide_free.fval(  :,1)); prob_high_glide_free.fval_CQ_sorted   = prob_high_glide_free.fval(  prob_high_glide_free.sort_index  ,2);
[prob_smooth_glide_free.fval_CP_sorted , prob_smooth_glide_free.sort_index] = sort(prob_smooth_glide_free.fval(:,1)); prob_smooth_glide_free.fval_CQ_sorted = prob_smooth_glide_free.fval(prob_smooth_glide_free.sort_index,2);
[prob_free_limit_case.fval_CP_sorted   , prob_free_limit_case.sort_index  ] = sort(prob_free_limit_case.fval(  :,1)); prob_free_limit_case.fval_CQ_sorted   = prob_free_limit_case.fval(  prob_free_limit_case.sort_index  ,2);


% Line colors
base_colors = lines(7);
our_colors = [0 , 0, 0, ; base_colors(1,:) ; base_colors(2,:); base_colors(5,:); base_colors(4,:)];
figure(1002);
h1 = plot(-prob_FFA_free.cf(1)                  , prob_FFA_free.cf(2)                  , 'o' , 'Color', our_colors( 1,:)); hold on; grid on; 
h2 = plot(-prob_FFA_free.fval_CP_sorted         , prob_FFA_free.fval_CQ_sorted         , '.-', 'Color', our_colors( 1,:));
h3 = plot(-prob_smooth_glide_free.fval_CP_sorted, prob_smooth_glide_free.fval_CQ_sorted, '.-', 'Color', our_colors( 2,:));
h4 = plot(-prob_high_glide_free.fval_CP_sorted  , prob_high_glide_free.fval_CQ_sorted  , '.-', 'Color', our_colors( 3,:));
h4b = plot(-prob_free_limit_case.fval_CP_sorted  , prob_free_limit_case.fval_CQ_sorted  , '.-', 'Color', our_colors(5,:));

legend('Original Design' , 'Feasible with FFA airfoil', 'Feasible with smooth glide airfoil', 'Feasible with high glide airfoil', 'Feasible with ideal airfoil', 'Location', 'NorthWest');
title('Optimization with probabilistic approach')
xlabel('C_P'); ylabel('CQ_{RBM}');
axis([0.41 0.51 0.45 0.70])

print -depsc fig/pareto_front_full_figure_no_isocosts.eps
savefig('fig/pareto_front_full_figure_no_isocosts.fig')


iso_cost_color_m = [0.4660    0.6740    0.1880]; %[0.1 0.4 0.1];
%iso_cost_color_p = [0.6740    0.4660    0.1880]; % base_colors(5,:); %[0.4 0.1 0.1];
%iso_cost_color_m = 0.4 * [1 1 1]; %[0.1 0.4 0.1];
iso_cost_color_p = 0.4 * [1 1 1]; % base_colors(5,:); %[0.4 0.1 0.1];
h5 = plot(CP_isocost_range, CQ_RBM_iso_cost_m10, '-.', 'Color', iso_cost_color_m ); hold on
plot(CP_isocost_range, CQ_RBM_iso_cost_m09, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m08, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m07, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m06, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m05, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m04, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m03, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m02, '-.', 'Color', iso_cost_color_m )
plot(CP_isocost_range, CQ_RBM_iso_cost_m01, '-.', 'Color', iso_cost_color_m )
h6 = plot(CP_isocost_range, CQ_RBM_iso_cost_000, '-.', 'Color', [0 0 0] );
h7 = plot(CP_isocost_range, CQ_RBM_iso_cost_p01, '-.', 'Color', iso_cost_color_p );
plot(CP_isocost_range, CQ_RBM_iso_cost_p02, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p03, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p04, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p05, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p06, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p07, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p08, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p09, '-.', 'Color', iso_cost_color_p )
plot(CP_isocost_range, CQ_RBM_iso_cost_p10, '-.', 'Color', iso_cost_color_p )
%plot(-stat_FFA_free.cf(1)         , stat_FFA_free.cf(2)         , 'ok'); % hold on; grid on; plot(-stat_FFA_free.fval_CP_sorted   , stat_FFA_free.fval_CQ_sorted   , '.k-'); 
%plot(-stat_high_glide_free.cf(1)  , stat_high_glide_free.cf(2)   , 'or'); plot(-stat_high_glide_free.fval_CP_sorted, stat_high_glide_free.fval_CQ_sorted, '.r-');
%axis([0.44 0.50 0.3 0.55])
%xlabel('Cp'); ylabel('CQ_{RBM}')
legend([h1 h2 h3 h4 h4b h5 h6 h7], 'Original Design' , ...
                               'Feasible with FFA airfoil', ...
                               'Feasible with smooth glide airfoil', ...
                               'Feasible with high glide airfoil', ...
                               'Feasible with ideal airfoil', ...
                               'Isocost (-0.4% LCOE steps)', ...
                               'Isocost (Innwind baseline)', ...
                               'Isocost (+0.4% LCOE steps)', ...
                               'Location', 'NorthWest');

print -depsc fig/pareto_front_full_figure_isocosts_0p4pc_v2.eps
savefig('fig/pareto_front_full_figure_isocosts_0p4pc_v2.fig')

set(gcf, 'PaperType', 'A5')
orient landscape 
print -dpdf fig/pareto_front_full_figure_isocosts_0p4pc_v2.dpf


% figure(101);
% plot(-stat_FFA_free.cf(1)         , stat_FFA_free.cf(2)         , 'ok'); hold on; grid on; plot(-stat_FFA_free.fval_CP_sorted   , stat_FFA_free.fval_CQ_sorted   , '.k-'); 
% plot(-stat_high_glide_free.cf(1)  , stat_high_glide_free.cf(2)  , 'or'); plot(-stat_high_glide_free.fval_CP_sorted, stat_high_glide_free.fval_CQ_sorted, '.r-');
% plot(-stat_smooth_glide_free.cf(1), stat_smooth_glide_free.cf(2), 'ob'); plot(-stat_smooth_glide_free.fval_CP_sorted, stat_smooth_glide_free.fval_CQ_sorted, '.b-');
% legend('FFA Design' , 'FFA Design', 'High glide', 'High glide', 'Smooth glide', 'Smooth glide'); title('Stat optimization')
% xlabel('C_P'); ylabel('CQ_{RBM');
% 




% 
% % Plot AEP vs cp
% cp_range = linspace(0.3, 16/27);
% plot(cp_range, AEP_low, cp_range, AEP_high, cp_range, AEP, 'k'); grid on





function [AEP, AEP_low, AEP_high] = AEP_internal_fun(ICM, CP)
    AEP_low  = ICM.a_low_cp  + ICM.b_low_cp  * CP;
    AEP_high = ICM.a_high_cp + ICM.b_high_cp * CP;
    
    AEP                         = zeros(size(CP));
    AEP(CP <  ICM.CP_threshold) = AEP_low( CP <  ICM.CP_threshold);
    AEP(CP >= ICM.CP_threshold) = AEP_high(CP >= ICM.CP_threshold);
end





