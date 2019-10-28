function [pareto_bh_sorted , pareto_ld_bar_sorted] = process_awec_ldtc_pareto_front(results_file)

results = load(results_file);

pareto_cf = results.results.fval;

% Extract Interesting Parameters from Pareto Front
pareto_bh = - pareto_cf(:,2);
pareto_ld_bar = 100 ./ pareto_cf(:,1);

% Now sort pareto front
[pareto_bh_sorted, pareto_bh_sorting_index] = sort(pareto_bh);
pareto_ld_bar_sorted         = pareto_ld_bar(        pareto_bh_sorting_index);

end