
function [cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper_custom_p(CPR, SPR, EC_cell, u_order, l_order)
% APE_FUNCTION_WRAPPER wraps the evaluation of a cell array of
% experimental_cases of airfoil_polar IDkind for a given simulation profile
% (SPR) by an airfoil_polar_evaluator (APE) object into a batchable job for
% distributed computing.

    % Create and airfoil polar evaluator
    APE = airfoil_polar_evaluator(CPR, SPR, EC_cell);
    % Set custom parametrization order
    APE.u_order = u_order;
    APE.l_order = l_order;
    % Initialize a local context
    APE.initialize_local_context();
    % Run all simulations
    APE.run_simulations();
    % Compare results (no penalty for non-convergence right now)
    [cl_global_accuracy , cm_global_accuracy , cd_global_accuracy] = APE.compare_simulation_results();
end