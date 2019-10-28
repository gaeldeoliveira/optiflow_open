function interpolated_values = polar_tensor_interpolator(polar_tensors, interpolated_tensor, alpha, Re, tc)
    % Sanitize inputs
    
    % Make all alpha values fit inside bounds
    alpha_safe_a = max(alpha       , min(polar_tensors.alpha_range+eps(10)));
    alpha_safe_b = min(alpha_safe_a, max(polar_tensors.alpha_range-eps(10)));
    % Make all Re values fit inside bounds
    Re_safe_a = max(Re             , min(polar_tensors.Re_range   +eps(10)));
    Re_safe_b = min(Re_safe_a      , max(polar_tensors.Re_range   -eps(10)));
    % Make all tc values fit inside bounds
    tc_safe_a = max(tc             , min(polar_tensors.tc_range   +eps(10)));
    tc_safe_b = min(tc_safe_a      , max(polar_tensors.tc_range   -eps(10)));
    % Interpolate
    interpolated_values = interpn(polar_tensors.alpha_tensor, polar_tensors.Re_tensor, polar_tensors.tc_tensor, interpolated_tensor, alpha_safe_b, Re_safe_b, tc_safe_b);
    
end