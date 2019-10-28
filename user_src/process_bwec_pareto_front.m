function [pareto_dCL_sorted , pareto_ld_bar_sorted] = process_bwec_pareto_front(results_file)

results = load(results_file);

pareto_cf = results.results.fval;

pareto_dCL                   =      - pareto_cf(:,1);
pareto_ld_bar                = 100 ./ pareto_cf(:,2);

% Now sort pareto front
[pareto_dCL_sorted, pareto_dCL_sorting_index] = sort(pareto_dCL);
pareto_ld_bar_sorted         = pareto_ld_bar(        pareto_dCL_sorting_index);

end