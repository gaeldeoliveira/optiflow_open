% A very simple BEM code bogo2_alfa7
addpath src
tic

% % Load moded Inwind rotor description
inwind_moded_RWT = load('../rotor_integration/graph_digitization/inwind_moded_RWT.mat'); RWT = inwind_moded_RWT.inwind_moded_RWT;
% % Load airfoil polar tensors
FFA_tensors = load('../rotor_integration/airfoil_families/FFA/FFA_free.mat'); polar_tensors = FFA_tensors.polar_tensors;
% % Set reference wind speed and rotational speed
U_ref = 8; lambda_ref = RWT.lambda_design;
Omega_ref = lambda_ref * U_ref/  RWT.R;
% % Define BemCase
% Scalar Inputs
BC.R       = RWT.R     ;     % [m     ]  Rotor tip  radius (corresponds to R)
BC.R_root  = RWT.R_hub ;     % [m     ]  Rotor root radius
BC.Omega   = Omega_ref ;     % [rad/s ]  Rotational spped
BC.U_inf   = U_ref     ;     % [m/s   ]  Inflow speed
BC.B       = 3         ;     % [#     ]  Number of blades
BC.nu_inf  = 1.46e-5   ;     % [m2/s  ]  Fluid Viscosity (nu = mu/rho)
BC.rho_inf = 1.225     ;     % [kg/m3 ]  Fluid Density
BC.N_bins  = 60        ;     % [#     ]  Number of computational elements along blade
% Functional Inputs
%BC.c_over_R_of_mu_fun         = @(mu ) 0.1 * ones(size(mu));                                       % Make a constant   chord distribution
BC.c_over_R_of_mu_fun         = @(mu ) interp1(RWT.r_chord/ RWT.R, RWT.val_chord / RWT.R, mu);      % Chord distribution
%BC.theta_deg_of_mu_fun        = @(mu ) atan2(BC.U_inf, 70 * (mu*BC.R)) * 180 /pi();% - 5 ;         % Make a helical    twist distribution (don't reinterpolate)
BC.theta_deg_of_mu_fun        = @(mu ) interp1(RWT.r_twist_deg/ RWT.R, RWT.val_twist_deg      , mu);      % Twist distribution
BC.tc_of_mu_fun               = @(mu ) interp1(RWT.r_thickness/ RWT.R, RWT.val_thickness / 100, mu);      % Thickness distribution
BC.cl_of_alpha_deg_and_re_fun = @(alpha, Re, tc) polar_tensor_interpolator(polar_tensors, polar_tensors.cl_tensor, alpha, Re, tc);
BC.cd_of_alpha_deg_and_re_fun = @(alpha, Re, tc) polar_tensor_interpolator(polar_tensors, polar_tensors.cd_tensor, alpha, Re, tc);

% Solution details
BC.bool_loss       = true;
BC.bool_glauert    = true;
BC.bool_tan_solve  = true;

BC.bool_prob_coef  = false;
BC.bool_exceedance = false;


[RES_axi, RES_tan, BD, BR] = bogo2_alfa7_BEMsolver(BC);

disp(['BR.P_rotor_annuli    = ' num2str(BR.P_rotor_annuli   )])
disp(['BR.P_rotor_segments  = ' num2str(BR.P_rotor_segments )])
disp(['BR.CP_rotor_annuli   = ' num2str(BR.CP_rotor_annuli  )])
disp(['BR.CP_rotor_segments = ' num2str(BR.CP_rotor_segments)])
disp(['BR.CT_rotor_annuli   = ' num2str(BR.CT_rotor_annuli  )])
disp(['BR.CT_rotor_segments = ' num2str(BR.CT_rotor_segments)])
disp(['BR.CQ_axi_RBM        = ' num2str(BR.CQ_axi_root_bending_moment)])
disp(['BR.CQ_tan_RBM        = ' num2str(BR.CQ_tan_root_bending_moment)])
disp(['BR.RES_axi           = ' num2str(max(BR.RES_axi)     )])
disp(['BR.RES_tan           = ' num2str(max(BR.RES_tan)     )])

% % Now make 
% [RES_axi ; RES_tan]
% BR.lambda
% CP_rotor_annuli
% CP_rotor_segments
toc


