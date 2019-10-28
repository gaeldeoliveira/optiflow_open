% A very simple BEM code bogo2_alfa7 hook to GAM optimization
% Driven from external driver script

% % Load moded Inwind rotor description
inwind_moded_RWT = load(PB.base_blade_geometry_file); RWT = inwind_moded_RWT.inwind_moded_RWT;
% % Load airfoil polar tensors
% Load airfoil polars
polar_tensors_container = load(PB.polar_tensors_file); polar_tensors = polar_tensors_container.polar_tensors;
% % Set reference wind speed and rotational speed
U_ref = 10.6; lambda_ref = RWT.lambda_design;
Omega_ref = lambda_ref * U_ref/  RWT.R;
% % Load angle of attack perturbation distribution
std_distribution = load(PB.std_distribution_file); std_distribution = std_distribution.std_distribution;

% % Define BemCase
% Scalar Inputs
BC.R       = RWT.R     ;     % [m     ]  Rotor tip  radius (corresponds to R)
BC.R_root  = RWT.R_hub ;     % [m     ]  Rotor root radius
BC.Omega   = Omega_ref ;     % [rad/s ]  Rotational spped
BC.U_inf   = U_ref     ;     % [m/s   ]  Inflow speed
BC.B       = 3         ;     % [#     ]  Number of blades
BC.nu_inf  = 1.46e-5   ;     % [m2/s  ]  Fluid Viscosity (nu = mu/rho)
BC.rho_inf = 1.225     ;     % [kg/m3 ]  Fluid Density
BC.N_bins  = 40        ;     % [#     ]  Number of computational elements along blade
% Functional Inputs (blade)
%BC.c_over_R_of_mu_fun         = @(mu ) 0.1 * ones(size(mu));                                                 % Make a constant   chord distribution
BC.c_over_R_of_mu_fun               = @(mu ) interp1(RWT.r_chord/ RWT.R, RWT.val_chord / RWT.R, mu);          % [adim.] Chord distribution (not in percent)
%BC.theta_deg_of_mu_fun        = @(mu ) atan2(BC.U_inf, 70 * (mu*BC.R)) * 180 /pi();% - 5 ;                   % Make a helical    twist distribution (don't reinterpolate)
BC.theta_deg_of_mu_fun              = @(mu ) interp1(RWT.r_twist_deg/ RWT.R, RWT.val_twist_deg    , mu);      % [deg  ] Twist distribution
BC.tc_of_mu_fun                     = @(mu ) interp1(RWT.r_thickness/ RWT.R, RWT.val_thickness/100, mu);      % [adim.] Thickness distribution (not in percent)
% Functional Inputs (aerodynamic coefficients)
BC.cl_of_alpha_deg_and_re_fun       = @(alpha, Re, tc) polar_tensor_interpolator(polar_tensors, polar_tensors.cl_tensor, alpha, Re, tc);
BC.cd_of_alpha_deg_and_re_fun       = @(alpha, Re, tc) polar_tensor_interpolator(polar_tensors, polar_tensors.cd_tensor, alpha, Re, tc);
% Functional Inputs (probabilistic aerodynamic coefficients)
BC.cl_of_alpha_deg_and_re_prob_fun  = @(alpha, Re, tc, mu) probabilistic_polar_tensor_interpolator(polar_tensors, polar_tensors.cl_tensor, std_distribution, alpha, Re, tc, mu);
BC.cd_of_alpha_deg_and_re_prob_fun  = @(alpha, Re, tc, mu) probabilistic_polar_tensor_interpolator(polar_tensors, polar_tensors.cd_tensor, std_distribution, alpha, Re, tc, mu);