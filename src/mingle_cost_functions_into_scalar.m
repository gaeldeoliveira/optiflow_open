function cf = mingle_cost_functions_into_scalar(cfg_val, mingle_factor)
    % Mingle multi-objective costfunction into scalar, deals with vectorized
    % inputs!
    
    % Extract single cost functions from global_cost_function values 
    cf1 = cfg_val(:,1);
    cf2 = cfg_val(:,2);
    % Mingle results
    cf  = cf1 .* (1 - mingle_factor) + mingle_factor .* cf2;
    % And return!
end