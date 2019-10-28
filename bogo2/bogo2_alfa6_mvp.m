% A very simple BEM code bogo2_alfa6
addpath src
tic
% % Load airfoil family

% % Define BemCase
% Scalar Inputs
BC.R       = 1         ;     % [m     ]  Rotor tip  radius (corresponds to R)
BC.R_root  = 0.2       ;     % [m     ]  Rotor root radius
BC.Omega   = 50      ;     % [rad/s ]  Rotational spped
BC.U_inf   = 10        ;     % [m/s   ]  Inflow speed
BC.B       = 2         ;     % [#     ]  Number of blades
BC.nu_inf  = 1.46e-5   ;     % [m2/s  ]  Fluid Viscosity (nu = mu/rho)
BC.rho_inf = 1.225     ;     % [kg/m3 ]  Fluid Density
BC.N_bins  = 10        ; 
% Functional Inputs
BC.c_over_R_of_mu_fun         = @(mu ) 0.1 * ones(size(mu));                                                 % Make a constant   chord distribution
%theta_deg_of_mu_fun       = @(mu ) 10  * ones(size(mu));                                                 % Make a constant   twist distribution
BC.theta_deg_of_mu_fun        = @(mu ) atan2(BC.U_inf, 70 * (mu*BC.R)) * 180 /pi() - 5 ;                        % Make a helical    twist distribution (don't reinterpolate)
BC.cl_of_alpha_deg_and_re_fun = @(alpha_deg, Re) 2*pi*sin(alpha_deg*pi/180).*cos(alpha_deg*pi/180);
%BC.cl_of_alpha_deg_and_re_fun = @(alpha_deg, Re) zeros(size(alpha_deg));
BC.cd_of_alpha_deg_and_re_fun = @(alpha_deg, Re) 0.01 * (1+BC.cl_of_alpha_deg_and_re_fun(alpha_deg, Re)).^2;
%BC.cd_of_alpha_deg_and_re_fun = @(alpha_deg, Re) zeros(size(alpha_deg));
% Solution details
BC.bool_loss       = true;
BC.bool_glauert    = true;
% BC.bool_axi_solve  = true;
BC.bool_tan_solve  = true;


[RES_axi, RES_tan, BD, BR] = bogo2_alfa6_BEMsolver(BC);

disp(['BR.CP_rotor_annuli   = ' num2str(BR.CP_rotor_annuli  )])
disp(['BR.CP_rotor_segments = ' num2str(BR.CP_rotor_segments)])
disp(['BR.RES_axi           = ' num2str(max(BR.RES_axi)     )])
disp(['BR.RES_tan           = ' num2str(max(BR.RES_tan)     )])

% % Now make 
% [RES_axi ; RES_tan]
% BR.lambda
% CP_rotor_annuli
% CP_rotor_segments
toc

