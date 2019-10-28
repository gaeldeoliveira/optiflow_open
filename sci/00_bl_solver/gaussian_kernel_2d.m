function u_theta = gaussian_kernel_2d(r, sigma)
% 2d gaussian kernel for speed
u_theta = 1./(2*pi()*r) .* (1 - exp(-(r/sigma).^2));

end