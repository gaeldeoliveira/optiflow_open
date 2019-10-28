function [pareto_cm_sorted , pareto_ld_bar_sorted] = process_awec_ldcm_pareto_front(results_file)

results = load(results_file);

pareto_cf = results.results.fval;

% Extract Interesting Parameters from Pareto Front
pareto_cm      =    - pareto_cf(:,2);
pareto_ld_bar = 1 ./ pareto_cf(:,1);

% Now sort pareto front
[pareto_cm_sorted, pareto_cm_sorting_index] = sort(pareto_cm);
pareto_ld_bar_sorted         = pareto_ld_bar(pareto_cm_sorting_index);

end