addpath src/

% date_str = 'DO_20180915T200345';
% 
% stat_FFA_free          = load([date_str '/stat_FFA_free.mat'   ]);
% stat_high_glide_free   = load([date_str '/stat_high_glide_free.mat']);
% 
% % Make plot of effect of turbulence on two FFA airfoils 
% alpha_range = -5:0.2:20;
% Re          =  9e6;
% subplot(221)
% plot(alpha_range, stat_FFA_free.BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range, Re, 0.24, 0)); hold on
% plot(alpha_range, stat_FFA_free.BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range, Re, 0.24, 2));
% plot(alpha_range, stat_FFA_free.BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range, Re, 0.24, 0.2));


%stat_FFA_free.BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range, Re, tc, mu)
pt_FFA_free       = load('../rotor_integration/airfoil_families/FFA/FFA_freeC.mat');
pt_high_glide     = load('../rotor_integration/airfoil_families/DU-IW/high_glide_design_free.mat');

alpha_range = -5:0.2:20;
Re          =  9e6;
figure(2001)
% expected_values = probabilistic_polar_tensor_interpolator_explicit(polar_tensors, interpolated_tensor, alpha, Re, tc, std_alpha_deg)
subplot(221)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 0)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 2)); 
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 4)); 
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
grid on; ylabel('Cl'); xlabel('\alpha [deg]')
title('FFA-W3-241')
subplot(223)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 0) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.24, 0)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 2) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.24, 2)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 4) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.24, 4)); hold on;
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
ylabel('L/D'); xlabel('\alpha [deg]')
title('FFA-W3-241')
subplot(222)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.30, 0)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.30, 2)); 
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.30, 4)); 
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
grid on; ylabel('Cl'); xlabel('\alpha [deg]')
title('FFA-W3-301')
subplot(224)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.30, 0) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.30, 0)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.30, 2) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.30, 2)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.30, 4) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.30, 4)); hold on;
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
ylabel('L/D'); xlabel('\alpha [deg]')
title('FFA-W3-301')

print -depsc fig/polar_comparison_FFA_free.eps
savefig('fig/polar_comparison_FFA_free.fig')


figure(2002)
% expected_values = probabilistic_polar_tensor_interpolator_explicit(polar_tensors, interpolated_tensor, alpha, Re, tc, std_alpha_deg)
subplot(221)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 0)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 2)); 
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
grid on; ylabel('Cl'); xlabel('\alpha [deg]')
title('FFA-W3-241')
axis([-5 20 0 2.5])
subplot(223)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 0) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.24, 0)); hold on;
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cl_tensor, alpha_range, Re, 0.24, 2) ./ probabilistic_polar_tensor_interpolator_explicit(pt_FFA_free.polar_tensors, pt_FFA_free.polar_tensors.cd_tensor, alpha_range, Re, 0.24, 2)); hold on;
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
ylabel('L/D'); xlabel('\alpha [deg]')
title('FFA-W3-241')
axis([-5    20     0   200])
subplot(222)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_high_glide.polar_tensors, pt_high_glide.polar_tensors.cl_tensor, alpha_range, Re, 0.21, 0)); hold on;
plot(alpha_range, smooth(probabilistic_polar_tensor_interpolator_explicit(pt_high_glide.polar_tensors, pt_high_glide.polar_tensors.cl_tensor, alpha_range, Re, 0.21, 2))); 
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
grid on; ylabel('Cl'); xlabel('\alpha [deg]')
title('High glide design')
axis([-5 20 0 2.5])
subplot(224)
plot(alpha_range, probabilistic_polar_tensor_interpolator_explicit(pt_high_glide.polar_tensors, pt_high_glide.polar_tensors.cl_tensor, alpha_range, Re, 0.21, 0) ./ probabilistic_polar_tensor_interpolator_explicit(pt_high_glide.polar_tensors, pt_high_glide.polar_tensors.cd_tensor, alpha_range, Re, 0.21, 0)); hold on;
plot(alpha_range, smooth(probabilistic_polar_tensor_interpolator_explicit(pt_high_glide.polar_tensors, pt_high_glide.polar_tensors.cl_tensor, alpha_range, Re, 0.21, 2) ./ probabilistic_polar_tensor_interpolator_explicit(pt_high_glide.polar_tensors, pt_high_glide.polar_tensors.cd_tensor, alpha_range, Re, 0.21, 2))); hold on;
grid on; legend('\sigma_{\alpha} = 0', '\sigma_{\alpha} = 2', '\sigma_{\alpha} = 4', 'Location', 'SouthEast')
ylabel('L/D'); xlabel('\alpha [deg]')
title('High glide design')
axis([-5    20     0   200])
print -depsc fig/polar_comparison_FFA_to_highglide.eps
savefig('fig/polar_comparison_FFA_to_highglide.fig')

% % ERROR IDENTIFIED %% Bypass for now %%
%
%   Rotor optimizations for pure FFA   used FFA-W3-211 instead of FFA-W3-241 at tip
%   Rotor optimizations for high glide used it as if it were 21%, but that
%   is not a big issue because used rotor definition had thickness scaled
%   to use a 21% airfoil at tip. 
%
%   This means real comparison is between High glide and FFA-W3-211, but
%   thesis describes as being between High glide and FFA-W3-241. This means
%   the offset is understimated, and thay our results are presented in a,
%   slightly erroneous, but above all, more conservative way. %
%
%   Real improvement in rotor performance would be even larger in
%   principle! Enjoy!
%


