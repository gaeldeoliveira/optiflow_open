% List cases


date_str = 'DO_20180914T222208';

stat_FFA_free          = load([date_str '/stat_FFA_free_tmp.mat'   ]);
stat_high_glide_free   = load([date_str '/stat_high_glide_free.mat']);
stat_smooth_glide_free = load([date_str '/stat_smooth_glide_free.mat']);

prob_FFA_free          = load([date_str '/prob_FFA_free.mat'   ]);
prob_high_glide_free   = load([date_str '/prob_high_glide_free.mat']);
prob_smooth_glide_free = load([date_str '/prob_high_glide_free.mat']);



% % Post stuff
% Reorder pareto fronts (stat)
[stat_FFA_free.fval_CP_sorted          , stat_FFA_free.sort_index         ] = sort(stat_FFA_free.fval(         :,1)); stat_FFA_free.fval_CQ_sorted          = stat_FFA_free.fval(         stat_FFA_free.sort_index         ,2);
[stat_high_glide_free.fval_CP_sorted   , stat_high_glide_free.sort_index  ] = sort(stat_high_glide_free.fval(  :,1)); stat_high_glide_free.fval_CQ_sorted   = stat_high_glide_free.fval(  stat_high_glide_free.sort_index  ,2);
[stat_smooth_glide_free.fval_CP_sorted , stat_smooth_glide_free.sort_index] = sort(stat_smooth_glide_free.fval(:,1)); stat_smooth_glide_free.fval_CQ_sorted = stat_smooth_glide_free.fval(stat_smooth_glide_free.sort_index,2);
% Reorder pareto fronts (prob)
[prob_FFA_free.fval_CP_sorted          , prob_FFA_free.sort_index         ] = sort(prob_FFA_free.fval(         :,1)); prob_FFA_free.fval_CQ_sorted          = prob_FFA_free.fval(         prob_FFA_free.sort_index         ,2);
[prob_high_glide_free.fval_CP_sorted   , prob_high_glide_free.sort_index  ] = sort(prob_high_glide_free.fval(  :,1)); prob_high_glide_free.fval_CQ_sorted   = prob_high_glide_free.fval(  prob_high_glide_free.sort_index  ,2);
[prob_smooth_glide_free.fval_CP_sorted , prob_smooth_glide_free.sort_index] = sort(prob_smooth_glide_free.fval(:,1)); prob_smooth_glide_free.fval_CQ_sorted = prob_smooth_glide_free.fval(prob_smooth_glide_free.sort_index,2);


% Plot pareto front (stat)
figure(101);
plot(-stat_FFA_free.cf(1)         , stat_FFA_free.cf(2)         , 'ok'); hold on; grid on; plot(-stat_FFA_free.fval_CP_sorted   , stat_FFA_free.fval_CQ_sorted   , '.k-'); 
plot(-stat_high_glide_free.cf(1)  , stat_high_glide_free.cf(2)  , 'or'); plot(-stat_high_glide_free.fval_CP_sorted, stat_high_glide_free.fval_CQ_sorted, '.r-');
plot(-stat_smooth_glide_free.cf(1), stat_smooth_glide_free.cf(2), 'ob'); plot(-stat_smooth_glide_free.fval_CP_sorted, stat_smooth_glide_free.fval_CQ_sorted, '.b-');
legend('FFA Design' , 'FFA Design', 'High glide', 'High glide', 'Smooth glide', 'Smooth glide'); title('Stat optimization')
%axis([0.45 0.5054 0.40 0.55])
% Plot pareto front (prob)
figure(102);
plot(-prob_FFA_free.cf(1)         , prob_FFA_free.cf(2)         , 'ok'); hold on; grid on; plot(-prob_FFA_free.fval_CP_sorted   , prob_FFA_free.fval_CQ_sorted   , '.k-');
plot(-prob_high_glide_free.cf(1)  , prob_high_glide_free.cf(2)  , 'or'); plot(-prob_high_glide_free.fval_CP_sorted, prob_high_glide_free.fval_CQ_sorted, '.r-');
plot(-prob_smooth_glide_free.cf(1), prob_smooth_glide_free.cf(2), 'ob'); plot(-prob_smooth_glide_free.fval_CP_sorted, prob_smooth_glide_free.fval_CQ_sorted, '.b-');
legend('FFA Design' , 'FFA Design', 'High glide', 'High glide', 'Smooth glide', 'Smooth glide');  title('Prob optimization')
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