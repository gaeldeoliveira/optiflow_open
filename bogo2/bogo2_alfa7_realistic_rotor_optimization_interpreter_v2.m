% List cases
addpath src/

date_str = 'DO_20180915T200345';

stat_FFA_free          = load([date_str '/stat_FFA_free.mat'   ]);
stat_high_glide_free   = load([date_str '/stat_high_glide_free.mat']);
stat_smooth_glide_free = load([date_str '/stat_smooth_glide_free.mat']);
stat_free_limit_case   = load([date_str '/stat_FFA_free_limit_case.mat']);

prob_FFA_free          = load([date_str '/prob_FFA_free.mat'   ]);
prob_high_glide_free   = load([date_str '/prob_high_glide_free.mat']);
prob_smooth_glide_free = load([date_str '/prob_smooth_glide_free.mat']);
prob_free_limit_case   = load([date_str '/prob_FFA_free_limit_case.mat']);



% % Post stuff
% Reorder pareto fronts (stat)
[stat_FFA_free.fval_CP_sorted          , stat_FFA_free.sort_index         ] = sort(stat_FFA_free.fval(         :,1)); stat_FFA_free.fval_CQ_sorted          = stat_FFA_free.fval(         stat_FFA_free.sort_index         ,2);
[stat_high_glide_free.fval_CP_sorted   , stat_high_glide_free.sort_index  ] = sort(stat_high_glide_free.fval(  :,1)); stat_high_glide_free.fval_CQ_sorted   = stat_high_glide_free.fval(  stat_high_glide_free.sort_index  ,2);
[stat_smooth_glide_free.fval_CP_sorted , stat_smooth_glide_free.sort_index] = sort(stat_smooth_glide_free.fval(:,1)); stat_smooth_glide_free.fval_CQ_sorted = stat_smooth_glide_free.fval(stat_smooth_glide_free.sort_index,2);
[stat_free_limit_case.fval_CP_sorted   , stat_free_limit_case.sort_index  ] = sort(stat_free_limit_case.fval(  :,1)); stat_free_limit_case.fval_CQ_sorted   = stat_free_limit_case.fval(  stat_free_limit_case.sort_index  ,2);

% Reorder pareto fronts (prob)
[prob_FFA_free.fval_CP_sorted          , prob_FFA_free.sort_index         ] = sort(prob_FFA_free.fval(         :,1)); prob_FFA_free.fval_CQ_sorted          = prob_FFA_free.fval(         prob_FFA_free.sort_index         ,2);
[prob_high_glide_free.fval_CP_sorted   , prob_high_glide_free.sort_index  ] = sort(prob_high_glide_free.fval(  :,1)); prob_high_glide_free.fval_CQ_sorted   = prob_high_glide_free.fval(  prob_high_glide_free.sort_index  ,2);
[prob_smooth_glide_free.fval_CP_sorted , prob_smooth_glide_free.sort_index] = sort(prob_smooth_glide_free.fval(:,1)); prob_smooth_glide_free.fval_CQ_sorted = prob_smooth_glide_free.fval(prob_smooth_glide_free.sort_index,2);
[prob_free_limit_case.fval_CP_sorted   , prob_free_limit_case.sort_index  ] = sort(prob_free_limit_case.fval(  :,1)); prob_free_limit_case.fval_CQ_sorted   = prob_free_limit_case.fval(  prob_free_limit_case.sort_index  ,2);

% Plot pareto front (stat)
base_colors = lines(7);
our_colors = [0 , 0, 0, ; base_colors(1,:) ; base_colors(2,:); base_colors(5,:); base_colors(4,:)];
figure(101);
plot(-stat_FFA_free.cf(1)                  , stat_FFA_free.cf(2)                    , 'o' , 'Color', our_colors(1,:)); hold on; grid on; 
plot(-stat_FFA_free.fval_CP_sorted         , stat_FFA_free.fval_CQ_sorted           , '.-', 'Color', our_colors(1,:)); 
plot(-stat_smooth_glide_free.fval_CP_sorted, stat_smooth_glide_free.fval_CQ_sorted  , '.-', 'Color', our_colors(3,:));
plot(-stat_high_glide_free.fval_CP_sorted  , stat_high_glide_free.fval_CQ_sorted    , '.-', 'Color', our_colors(2,:));
plot(-stat_free_limit_case.fval_CP_sorted  , stat_free_limit_case.fval_CQ_sorted    , '.-', 'Color', our_colors(4,:));
legend('Original Design' , 'Feasible with FFA airfoil', 'Feasible with smooth glide airfoil', 'Feasible with high glide airfoil', 'Feasible with ideal airfoil', 'Location', 'NorthWest');
title('Optimization with static approach')
xlabel('C_P'); ylabel('CQ_{RBM}');
axis([0.42 0.52 0.30 0.55])

print -depsc fig/pareto_front_full_figure_no_isocosts_stat_from_interpreter.eps
savefig('fig/pareto_front_full_figure_no_isocosts_stat_from_interpreter.fig')
set(gcf, 'PaperType', 'A5')
orient landscape
print -dpdf fig/pareto_front_full_figure_no_isocosts_stat_from_interpreter.pdf

%axis([0.45 0.5054 0.40 0.55])
% Plot pareto front (prob)
figure(102);
plot(-prob_FFA_free.cf(1)                  , prob_FFA_free.cf(2)                  , 'o' , 'Color', our_colors(1,:)); hold on; grid on; 
plot(-prob_FFA_free.fval_CP_sorted         , prob_FFA_free.fval_CQ_sorted         , '.-', 'Color', our_colors(1,:)); 
plot(-prob_smooth_glide_free.fval_CP_sorted, prob_smooth_glide_free.fval_CQ_sorted, '.-', 'Color', our_colors(3,:));
plot(-prob_high_glide_free.fval_CP_sorted  , prob_high_glide_free.fval_CQ_sorted  , '.-', 'Color', our_colors(2,:));
plot(-prob_free_limit_case.fval_CP_sorted  , prob_free_limit_case.fval_CQ_sorted  , '.-', 'Color', our_colors(4,:));
legend('Original Design' , 'Feasible with FFA airfoil', 'Feasible with smooth glide airfoil', 'Feasible with high glide airfoil', 'Feasible with ideal airfoil', 'Location', 'NorthWest');
title('Optimization with probabilistic approach')
xlabel('C_P'); ylabel('CQ_{RBM}');
axis([0.42 0.52 0.45 0.70])

print -depsc fig/pareto_front_full_figure_no_isocosts_prob_from_interpreter.eps
savefig('fig/pareto_front_full_figure_no_isocosts_prob_from_interpreter.fig')
set(gcf, 'PaperType', 'A5')
orient landscape
print -dpdf fig/pareto_front_full_figure_no_isocosts_prob_from_interpreter.pdf


base_colors = lines(7);
our_colors = [0 , 0, 0, ; base_colors(1,:) ; base_colors(2,:); base_colors(5,:); base_colors(4,:)];
subplot(221);
plot(-stat_FFA_free.cf(1)                  , stat_FFA_free.cf(2)                    , 'o' , 'Color', our_colors(1,:)); hold on; grid on; 
plot(-stat_FFA_free.fval_CP_sorted         , stat_FFA_free.fval_CQ_sorted           , '.-', 'Color', our_colors(1,:)); 
plot(-stat_smooth_glide_free.fval_CP_sorted, stat_smooth_glide_free.fval_CQ_sorted  , '.-', 'Color', our_colors(3,:));
plot(-stat_high_glide_free.fval_CP_sorted  , stat_high_glide_free.fval_CQ_sorted    , '.-', 'Color', our_colors(2,:));
plot(-stat_free_limit_case.fval_CP_sorted  , stat_free_limit_case.fval_CQ_sorted    , '.-', 'Color', our_colors(4,:));
legend('Original Design' , 'Feasible with FFA airfoil', 'Feasible with smooth glide airfoil', 'Feasible with high glide airfoil', 'Feasible with ideal airfoil', 'Location', 'NorthWest');
title('Optimization with static approach')
xlabel('C_P'); ylabel('CQ_{RBM}');
axis([0.42 0.52 0.30 0.55])


%axis([0.45 0.5054 0.40 0.55])
% Plot pareto front (prob)
subplot(222);
plot(-prob_FFA_free.cf(1)                  , prob_FFA_free.cf(2)                  , 'o' , 'Color', our_colors(1,:)); hold on; grid on; 
plot(-prob_FFA_free.fval_CP_sorted         , prob_FFA_free.fval_CQ_sorted         , '.-', 'Color', our_colors(1,:)); 
plot(-prob_smooth_glide_free.fval_CP_sorted, prob_smooth_glide_free.fval_CQ_sorted, '.-', 'Color', our_colors(3,:));
plot(-prob_high_glide_free.fval_CP_sorted  , prob_high_glide_free.fval_CQ_sorted  , '.-', 'Color', our_colors(2,:));
plot(-prob_free_limit_case.fval_CP_sorted  , prob_free_limit_case.fval_CQ_sorted  , '.-', 'Color', our_colors(4,:));
legend('Original Design' , 'Feasible with FFA airfoil', 'Feasible with smooth glide airfoil', 'Feasible with high glide airfoil', 'Feasible with ideal airfoil', 'Location', 'NorthWest');
title('Optimization with probabilistic approach')
xlabel('C_P'); ylabel('CQ_{RBM}');
axis([0.42 0.52 0.45 0.70])



% Manual positioning of legend here
print -depsc fig/pareto_front_composite_figure_no_isocosts_from_interpreter.eps
savefig('fig/pareto_front_composite_figure_no_isocosts_from_interpreter.fig')






%axis([0.45 0.5054 0.40 0.55])


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % Plot effect aoa perturbations on polar curve of expected values
% alpha_range = -5:0.2:20;
% figure;
% subplot(2,2,1)
% plot(alpha_range, BC.cl_of_alpha_deg_and_re_fun(     alpha_range,9e6,0.21  )); hold on; xlabel('\alpha [deg]')
% plot(alpha_range, BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range,9e6,0.21,1)); grid on; ylabel('C_l')
% subplot(2,2,2)
% plot(alpha_range, BC.cd_of_alpha_deg_and_re_fun(     alpha_range,9e6,0.21  )); hold on; xlabel('\alpha [deg]')
% plot(alpha_range, BC.cd_of_alpha_deg_and_re_prob_fun(alpha_range,9e6,0.21,1)); grid on; ylabel('C_d')
% subplot(2,2,3)
% plot(alpha_range, BC.cl_of_alpha_deg_and_re_fun(     alpha_range,9e6,0.21  ) ./ BC.cd_of_alpha_deg_and_re_fun(     alpha_range,3e6,0.21  )); hold on; xlabel('\alpha [deg]')
% plot(alpha_range, BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range,9e6,0.21,1) ./ BC.cd_of_alpha_deg_and_re_prob_fun(alpha_range,3e6,0.21,1)); grid on; ylabel('L/D')