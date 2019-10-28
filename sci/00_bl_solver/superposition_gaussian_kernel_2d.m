function u_theta = superposition_gaussian_kernel_2d(y_target, y_origin_vector, gamma_vector, sigma_vector)
% Identify number of sources
n_sources = length(y_origin_vector);
% Initialize output
u_theta = 0;

for n_source = 1:n_sources
    r_current        = y_target - y_origin_vector(n_source);
    gamma_current    = gamma_vector(n_source);
    sigma_current    = sigma_vector(n_source);
    kernel_current   = gaussian_kernel_2d(r_current, sigma_current);
    u_theta_current  = gamma_current * kernel_current;
    u_theta          = u_theta + u_theta_current;
end

end