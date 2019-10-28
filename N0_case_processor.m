% Add paths
addpath src/
addpath dev/
addpath user_src/

bool_recompute_data_file = false;
postprocessing_data_file = 'probabilistic/new_results/postprocessing_data_file';
postprocessing_N_cores   = 1;
postprocessing_date_str  = '20180910T013303'; % Only used for loading


if bool_recompute_data_file == true
    % Load sigma 20 case
    case_S20 = load('probabilistic/new_results/N0_new_prob_case_finalrun_sigma20_fullcase_20180909T233741.mat');
    % Migrate it
    case_S20.SC.N_cores = postprocessing_N_cores; 
    case_S20.SC.migrate_context('a');
    % Sort results structure
    case_S20.results_sorted = sort_pareto_front(case_S20.results);
    
    % Load sigma 00 case
    case_S00 = load('probabilistic/new_results/N0_new_prob_case_finalrun_sigma00_fullcase_20180909T173835.mat');
    % Migrate it
    case_S00.SC.N_cores = postprocessing_N_cores;
    case_S00.SC.migrate_context('b');
    % Sort results structure
    case_S00.results_sorted = sort_pareto_front(case_S00.results);
    
    % Now compute pareto fronts at sigma 20
    fval_S20_at_S20 = case_S20.GM.call_cost_function(case_S20.results_sorted.x); experiment_results_S20_at_S20 = case_S20.GM.CFG.last_experiment_set;
    fval_S20_at_S00 = case_S00.GM.call_cost_function(case_S20.results_sorted.x); experiment_results_S20_at_S00 = case_S00.GM.CFG.last_experiment_set;
    % And compute pareto fronts at sigma 00
    fval_S00_at_S20 = case_S20.GM.call_cost_function(case_S00.results_sorted.x); experiment_results_S00_at_S20 = case_S20.GM.CFG.last_experiment_set;
    fval_S00_at_S00 = case_S00.GM.call_cost_function(case_S00.results_sorted.x); experiment_results_S00_at_S00 = case_S00.GM.CFG.last_experiment_set;

    % Save post-processed data to   file
    save([postprocessing_data_file , '_' , datestr(now, 30), '.mat'])
else
    % Load post-processed data from file
    load([postprocessing_data_file , '_' , postprocessing_date_str, '.mat']);
end

% Filter out spurious points from S20 pareto front
S20_ok_index = not(fval_S20_at_S20(:,2) == 00);

% Now plot pareto fronts side by side
subplot(2,2,2)
plot(-100 * [fval_S00_at_S20(1,1) ; fval_S00_at_S20(:,1)], -100 * [0 ; fval_S00_at_S20(:,2)]); hold on
plot(-100 * fval_S20_at_S20(S20_ok_index,1), -100 * fval_S20_at_S20(S20_ok_index,2)); grid on

legend('S20 at S20', 'S00 at S20')
xlabel('max ((L^{bar}/D^{bar})_{free})'); ylabel('max ((L^{bar}/D^{bar})_{trip})'); 
legend('Optimized for \sigma = 0deg', 'Optimized for \sigma = 2deg', 'Location', 'SouthWest')
title('Feasibility frontier with \sigma=2deg')
axis([150 300 40 120])

subplot(2,2,1)
plot(-100 * [fval_S00_at_S00(1,1) ; fval_S00_at_S00(:           ,1) ], -100 * [0 ; fval_S00_at_S00(:,2)]); hold on
plot(-100 * [fval_S20_at_S00(1,1) ; fval_S20_at_S00(S20_ok_index,1) ], -100 * [0 ; fval_S20_at_S00(S20_ok_index,2)]); grid on

legend('S20 at S00', 'S00 at S00')
xlabel('max ((L^{bar}/D^{bar})_{free})'); ylabel('max ((L^{bar}/D^{bar})_{trip})'); 
legend('Optimized for \sigma = 0deg', 'Optimized for \sigma = 2deg', 'Location', 'SouthWest')
title('Feasibility frontier with \sigma=0deg')
axis([150 300 40 120])

% Chosen airfoils
n_DUWP2_S20A = 11; % 22;
n_DUWP2_S20B = 26;
n_DUWP2_S00A = 17;

coord_DUWP2_S20A =  [experiment_results_S20_at_S20{n_DUWP2_S20A,1}.coordinates.tx,experiment_results_S20_at_S20{n_DUWP2_S20A,1}.coordinates.tz];
coord_DUWP2_S20B =  [experiment_results_S20_at_S20{n_DUWP2_S20B,1}.coordinates.tx,experiment_results_S20_at_S20{n_DUWP2_S20B,1}.coordinates.tz];
coord_DUWP2_S00A =  [experiment_results_S00_at_S20{n_DUWP2_S00A,1}.coordinates.tx,experiment_results_S00_at_S20{n_DUWP2_S00A,1}.coordinates.tz];

save('probabilistic/new_results/coord_DUWP2_S20A2.air', 'coord_DUWP2_S20A', '-ascii', '-double');
save('probabilistic/new_results/coord_DUWP2_S20B2.air', 'coord_DUWP2_S20B', '-ascii', '-double');
save('probabilistic/new_results/coord_DUWP2_S00A2.air', 'coord_DUWP2_S00A', '-ascii', '-double');

% Now plot polars of these airfoils
subplot(2,2,1)
alpha_range = -1:0.2:14;
plot(alpha_range, experiment_results_S00_at_S20{n_DUWP2_S00A, 1}.ap.field_alpha_prob('cl', alpha_range, 0)); hold on
plot(alpha_range, experiment_results_S00_at_S20{n_DUWP2_S00A, 1}.ap.field_alpha_prob('cl', alpha_range, 2)); grid on
legend('Free transition, \sigma = 0deg', 'Free transition, \sigma = 2deg', 'Location', 'SouthEast')
axis([0 12 1 2.6]); yticks([1 1.4 1.8 2.2 2.6]);
title('DU-240-WP2-S00')

subplot(2,2,2)
alpha_range = -1:0.2:14;
plot(alpha_range, experiment_results_S20_at_S20{n_DUWP2_S20A, 1}.ap.field_alpha_prob('cl', alpha_range, 0)); hold on
plot(alpha_range, experiment_results_S20_at_S20{n_DUWP2_S20A, 1}.ap.field_alpha_prob('cl', alpha_range, 2)); grid on
legend('Free transition, \sigma = 0deg', 'Free transition, \sigma = 2deg', 'Location', 'SouthEast')
axis([0 12 1 2.6]); yticks([1 1.4 1.8 2.2 2.6]);
title('DU-240-WP2-S20')

subplot(2,2,3)
alpha_range = -1:0.2:14;
plot(alpha_range, experiment_results_S00_at_S20{n_DUWP2_S00A, 1}.ap.field_alpha_prob('cl', alpha_range, 0)./experiment_results_S00_at_S20{n_DUWP2_S00A, 1}.ap.field_alpha_prob('cd', alpha_range, 0)); hold on
plot(alpha_range, experiment_results_S00_at_S20{n_DUWP2_S00A, 1}.ap.field_alpha_prob('cl', alpha_range, 2)./experiment_results_S00_at_S20{n_DUWP2_S00A, 1}.ap.field_alpha_prob('cd', alpha_range, 2)); grid on
legend('Free transition, \sigma = 0deg', 'Free transition, \sigma = 2deg', 'Location', 'SouthEast')
axis([0 12 0 300]); % yticks([1 1.4 1.8 2.2 2.6]);
title('DU-240-WP2-S00')

subplot(2,2,4)
alpha_range = -1:0.2:14;
plot(alpha_range, experiment_results_S20_at_S20{n_DUWP2_S20A, 1}.ap.field_alpha_prob('cl', alpha_range, 0)./experiment_results_S20_at_S20{n_DUWP2_S20A, 1}.ap.field_alpha_prob('cd', alpha_range, 0)); hold on
plot(alpha_range, experiment_results_S20_at_S20{n_DUWP2_S20A, 1}.ap.field_alpha_prob('cl', alpha_range, 2)./experiment_results_S20_at_S20{n_DUWP2_S20A, 1}.ap.field_alpha_prob('cd', alpha_range, 2)); grid on
legend('Free transition, \sigma = 0deg', 'Free transition, \sigma = 2deg', 'Location', 'SouthEast')
axis([0 12 0 300]); % yticks([1 1.4 1.8 2.2 2.6]);
title('DU-240-WP2-S20')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%                                                                         %
%       % Prototype new code %                                            %
%                                                                         %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% % % Prototype new code
% ap = experiment_results_S00_at_S20{n_DUWP2_S00A, 1}.ap;
% 
% field_name    = 'cl';
% cl_data    = ap.raw_data.(field_name);
% alpha_data    = ap.raw_data.alpha;
% cl_safe    = cl_data;
% alpha_safe    = alpha_data;
% n_safe        = 0;
% 
% for n_data = 1:length(cl_data)
%     n_safe = n_safe + 1;
%     % Only introduce forward datapoint if cl>cl_max
%     if alpha_data(n_data) > ap.alpha_cl_max_local_limited
%         if (alpha_data(n_data) - alpha_data(n_data-1)) > 0.6 + eps(10)
%           alpha_safe = [alpha_safe(1:(n_safe-1)) ; alpha_data(n_data-1) + 0.6; alpha_safe(n_safe:end)];
%           cl_safe = [cl_safe(1:(n_safe-1)) ; cl_data(n_data)      ; cl_safe(n_safe:end)];
%           n_safe = n_safe + 1;
%         end
%     end
% end
% 
% cl_data    = ap.raw_data.cl;
% cd_data    = ap.raw_data.cd;
% alpha_data = ap.raw_data.alpha;
% cl_safe    = cl_data;
% cd_safe    = cd_data;
% alpha_safe = alpha_data;
% n_safe     = 0;
% 
% for n_data = 1:length(alpha_data)
%     n_safe = n_safe + 1;
%     % Only introduce forward datapoint if cl>cl_max
%     if alpha_data(n_data) > ap.alpha_cl_max_local_limited
%         if (alpha_data(n_data) - alpha_data(n_data-1)) > 0.6 + eps(10)
%           alpha_safe = [alpha_safe(1:(n_safe-1)) ; alpha_data(n_data-1) + 0.6; alpha_safe(n_safe:end)];
%           cl_safe    = [cl_safe(1:(n_safe-1)) ; cl_data(n_data)      ; cl_safe(n_safe:end)];
%           cd_safe    = [cd_safe(1:(n_safe-1)) ; cd_data(n_data)      ; cd_safe(n_safe:end)];
%           n_safe = n_safe + 1;
%         end
%     end
% end
% 
% subplot(212)
% plot(alpha_data, cl_data, 'x-'); hold on
% plot(alpha_safe, cl_safe, 'o-');
% subplot(211)
% plot(cd_data   , cl_data, 'x-'); hold on
% plot(cd_safe   , cl_safe, 'o-');




function results_sorted = sort_pareto_front(results)
    % Copy results structure
    results_sorted            = results;
    
    % Find sorting index based on first column of fval
    [~ , sort_index]          = sort(results.fval(:,1));
    
    % Sort fval in output structure
    results_sorted.fval       = results.fval(sort_index, :);
    
    % Sort x    in output structure
    results_sorted.x          = results.x(    sort_index, :);
    
    % Append sorting index for later reference
    results_sorted.sort_index = sort_index;

end