
stat_result = load('bogo2_alfa7_optim_20180907T023959.mat'     );
prob_result = load('bogo2_alfa7_prob_optim_20180907T023604.mat');


% Get pareto fronts and reorder them
fval_stat   = stat_result.fval; [~ , i_stat] = sort(fval_stat(:,1)); fval_stat = fval_stat(i_stat, :);
fval_prob   = prob_result.fval; [~ , i_prob] = sort(fval_prob(:,1)); fval_prob = fval_prob(i_prob, :);

% Make a colormap
clm = lines(2); blue = clm(1,:); red = clm(2,:);
% And plot !
plot(prob_result.cf(1), prob_result.cf(2), 'x', 'Color', blue);   hold on;  plot(fval_prob(:,1), fval_prob(:,2), 'x-', 'Color', blue);
plot(stat_result.cf(1), stat_result.cf(2), '+', 'Color', red );   hold on;  plot(fval_stat(:,1), fval_stat(:,2), '+-', 'Color', red );
grid on; legend('Stat: reference', 'Stat: Pareto', 'Prob: reference', 'Prob: Pareto'); xlabel('CP'); ylabel('CQ_root_bending_moment');





