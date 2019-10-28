function [cf1_rinv_sorted , cf2_rinv_sorted, x_sorted, cf2_sort_index] = reorder_and_inv_pareto_front(results, inv_fun_cf1, inv_fun_cf2)
    % Gets pareto front from results structure, sorts it in phenotype space, 
    % inverts post function and gathers corresponding genotypes into sorted
    % array.
    %
    % Inputs:
    %   results
    %   inv_fun_cf1     - handle to inverse of 1st goal CF.post_function
    %   inv_fun_cf2     - handle to inverse of 2nd goal CF.post_function
    %
    % Outputs:
    %   cf1_rinv_sorted - sorted array of pareto front phenotypes (1st cost function, sorted by decreasing cf2, and inverted with inv_fun_cf1)
    %   cf2_rinv_sorted - sorted array of pareto front phenotypes (2nd cost function, sorted by decreasing cf2, and inverted with inv_fun_cf2)
    %   x_sorted        - sorted array of pareto front entries genotypes
    %   cf2_sort_index  - sort index for correspondence to original results structure entries (x for genotypes, fval for phenotypes)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cf1_post = results.fval(:,1);
    cf2_post = results.fval(:,2);
    % Sort
    [~, cf2_sort_index] = sort(-cf2_post);
    % Sort post phenotypes
    cf1_post_sorted = cf1_post(cf2_sort_index);
    cf2_post_sorted = cf2_post(cf2_sort_index);
    % Reinverse post phenotypes
    cf1_rinv_sorted = inv_fun_cf1(cf1_post_sorted);
    cf2_rinv_sorted = inv_fun_cf2(cf2_post_sorted);
    % Gather genotypes
    x_sorted        = results.x(cf2_sort_index,:);
    % Return!
end