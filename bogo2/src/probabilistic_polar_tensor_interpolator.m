function expected_values = probabilistic_polar_tensor_interpolator(polar_tensors, interpolated_tensor, std_distribution, alpha, Re, tc, mu)
    % Sanitize inputs
    % Make all alpha values fit inside bounds
    mu_safe_a = max(mu       , min(std_distribution.mu_range+eps(1)));
    mu_safe_b = min(mu_safe_a, max(std_distribution.mu_range-eps(1)));
    
    % Find standard deviations
    std_alpha_deg = interp1(std_distribution.mu_range, std_distribution.std_alpha_deg_total, mu_safe_b);
    % Make relative angle of attack offsets
    alpha_offset_m4    = alpha - 2.00 * std_alpha_deg;
    alpha_offset_m3    = alpha - 1.50 * std_alpha_deg;
    alpha_offset_m2    = alpha - 1.00 * std_alpha_deg;
    alpha_offset_m1    = alpha - 0.50 * std_alpha_deg;
    alpha_offset_00    = alpha                       ;
    alpha_offset_p1    = alpha + 0.50 * std_alpha_deg;
    alpha_offset_p2    = alpha + 1.00 * std_alpha_deg;
    alpha_offset_p3    = alpha + 1.50 * std_alpha_deg;
    alpha_offset_p4    = alpha + 2.00 * std_alpha_deg;
    % Make angle of attack weight offsets
    weight_m4 = 0.0278; % cdf('norm', -1.75, 0, 1) - cdf('norm', -2.25, 0, 1);
    weight_m3 = 0.0656; % cdf('norm', -1.25, 0, 1) - cdf('norm', -1.75, 0, 1);
    weight_m2 = 0.1210; % cdf('norm', -0.75, 0, 1) - cdf('norm', -1.25, 0, 1);
    weight_m1 = 0.1747; % cdf('norm', -0.25, 0, 1) - cdf('norm', -0.75, 0, 1);
    weight_00 = 0.1974; % cdf('norm',  0.25, 0, 1) - cdf('norm', -0.25, 0, 1);
    weight_p1 = 0.1747; % cdf('norm',  0.75, 0, 1) - cdf('norm',  0.25, 0, 1);
    weight_p2 = 0.1210; % cdf('norm',  1.25, 0, 1) - cdf('norm',  0.75, 0, 1);
    weight_p3 = 0.0656; % cdf('norm',  1.75, 0, 1) - cdf('norm',  1.25, 0, 1);
    weight_p4 = 0.0278; % cdf('norm',  2.25, 0, 1) - cdf('norm',  1.75, 0, 1);
    % Compute total weights (to correct systematic bias of numerical
    % integration on finite bounds)
    total_weight = weight_m4 + weight_m3 + weight_m2 + weight_m1 + weight_00 + weight_p1 + weight_p2 + weight_p3 + weight_p4;
    % Evaluate expected interpolated tensor at all offsets
    interpolated_values_m4 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_m4, Re, tc);
    interpolated_values_m3 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_m3, Re, tc);
    interpolated_values_m2 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_m2, Re, tc);
    interpolated_values_m1 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_m1, Re, tc);
    interpolated_values_00 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_00, Re, tc);
    interpolated_values_p1 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_p1, Re, tc);
    interpolated_values_p2 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_p2, Re, tc);
    interpolated_values_p3 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_p3, Re, tc);
    interpolated_values_p4 = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha_offset_p4, Re, tc);
    % Evaluate expected value of coefficients
    expected_values        = ( weight_m4 * interpolated_values_m4 + ...
                               weight_m3 * interpolated_values_m3 + ...
                               weight_m2 * interpolated_values_m2 + ...
                               weight_m1 * interpolated_values_m1 + ...
                               weight_00 * interpolated_values_00 + ...
                               weight_p1 * interpolated_values_p1 + ...
                               weight_p2 * interpolated_values_p2 + ...
                               weight_p3 * interpolated_values_p3 + ...
                               weight_p4 * interpolated_values_p4 ) / total_weight;                
    
    % Done!
end