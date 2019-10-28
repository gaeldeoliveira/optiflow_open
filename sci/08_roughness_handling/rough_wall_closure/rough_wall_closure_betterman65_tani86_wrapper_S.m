function y = rough_wall_closure_betterman65_tani86_wrapper_S(x, Re_h, lambda)
        % Non-vectorized wrapper for finding (cf_input and Pi_wake) with
        % (H12 and Re_theta_S) as fixed inputs (targets) 
         
        % Demultiplex inputs
        cf_input = x(1);
        Pi_wake  = x(2);
        
        % Compute
        [H12, Re_theta_S, Re_theta_R] = rough_wall_closure_betterman65_tani86(cf_input, Pi_wake, Re_h, lambda);
        
        % Multiplex outputs
        y(1)     = H12       ;
        y(2)     = Re_theta_S / 1000;
end