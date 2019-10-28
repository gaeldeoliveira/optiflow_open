% A very simple BEM code bogo2_alfa1

% % Define BemCase
% Scalar Inputs
BC.R       = 1         ;     % [m     ]  Rotor tip  radius (corresponds to R)
BC.R_root  = 0.2       ;     % [m     ]  Rotor root radius
BC.Omega   = 1000        ;     % [rad/s ]  Rotational spped
BC.U_inf   = 10        ;     % [m/s   ]  Inflow speed
BC.B       = 2         ;     % [#     ]  Number of blades
BC.nu_inf  = 1.46e-5   ;     % [m2/s  ]  Fluid Viscosity (nu = mu/rho)
BC.rho_inf = 1.225     ;     % [kg/m3 ]  Fluid Density
BC.N_bins  = 10        ; 

BC.bool_loss       = true;
BC.bool_glauert    = true;
BC.bool_axi_solve  = true;
BC.bool_tan_solve  = true;

% Functional Inputs
BC.c_over_R_of_mu_fun         = @(mu ) 0.1 * ones(size(mu));                                                 % Make a constant   chord distribution
%theta_deg_of_mu_fun       = @(mu ) 10  * ones(size(mu));                                                 % Make a constant   twist distribution
BC.theta_deg_of_mu_fun        = @(mu ) atan2(U_inf, 70 * (mu*R)) * 180 /pi();%- 5 ;                        % Make a helical    twist distribution (don't reinterpolate)
BC.cl_of_alpha_deg_and_re_fun = @(alpha_deg, Re) 2*pi*sin(alpha_deg*pi/180).*cos(alpha_deg*pi/180);
BC.cd_of_alpha_deg_and_re_fun = @(alpha_deg, Re) 0.01 * (1+cl_of_alpha_deg_and_re_fun(alpha_deg, Re)).^2;
BC.cd_of_alpha_deg_and_re_fun = @(alpha_deg, Re) zeros(size(alpha_deg));


% % Parse BemCase fields

%

% Inputs
lambda = Omega * R / U_inf;

% Discretization
mu_bin_edges_vector = linspace(R_root/R,1,N_bins);                                                        % [adim ] Edges of each discretization element
mu_vector           = 0.5 * (mu_bin_edges_vector(1:(end-1)) + mu_bin_edges_vector(2:end));                % [adim ] Centers of each discretization element
bin_width_vector    = R    * (mu_bin_edges_vector(2:end) - mu_bin_edges_vector(1:(end-1)));               % [m    ] Width of each discretization element (blade segment)
bin_area_vector     = pi() * R.^2 * (mu_bin_edges_vector(2:end).^2 - mu_bin_edges_vector(1:(end-1)).^2);  % [m2   ] Area  of each discretization element (rotor annuli )

% Primary Variables (first guess/input)
a_axi = zeros(size(mu_vector));
a_tan = zeros(size(mu_vector));

for n = 1:100000
% Secondary variables: Positions and Sizes
r_vector         = mu_vector * R;
c_vector         = c_over_R_of_mu_fun(mu_vector) * R;
theta_deg        = theta_deg_of_mu_fun(mu_vector);

% Secondary variables: Speeds
U_axi = U_inf * (1 - a_axi);
U_tan = Omega * r_vector .* (1 + a_tan);
W     = sqrt(U_axi.^2 + U_tan.^2);

% Secondary variables: Reynolds number
Re_vector = W .* c_vector / nu_inf;

% Secondary variables: Angles
phi_rad   = atan2(U_axi, U_tan);
phi_deg   = phi_rad * 180 /pi();

alpha_deg = phi_deg - theta_deg;

% Secondary variables: Airfoil coefficients
cl_vector            = cl_of_alpha_deg_and_re_fun(alpha_deg, Re_vector);
cd_vector            = cd_of_alpha_deg_and_re_fun(alpha_deg, Re_vector);

% Secondary variables: Airfoil forces
L_per_unit_span      = 0.5 * rho_inf * W.^2 .* c_vector .* cl_vector;
D_per_unit_span      = 0.5 * rho_inf * W.^2 .* c_vector .* cd_vector;

% Secondary variables: Blade section forces
F_axi_per_unit_span  =   L_per_unit_span .* cos(phi_rad) + D_per_unit_span .* sin(phi_rad);
F_tan_per_unit_span  = - L_per_unit_span .* sin(phi_rad) + D_per_unit_span .* cos(phi_rad);

% Secondary variables: Blade bin forces
F_axi_blade_segments = F_axi_per_unit_span .* bin_width_vector;
F_tan_blade_segments = F_tan_per_unit_span .* bin_width_vector;

% Secondary variables: Rotor bin forces
F_axi_rotor_segments = F_axi_blade_segments * B;
F_tan_rotor_segments = F_tan_blade_segments * B;

% Secondary variables: Rotor bin moment
Q_tan_rotor_segments = F_tan_rotor_segments .* r_vector;

% Secondary variables: Tip loss factor
f_tip  = (B/2   ) * (R - r_vector   ) ./ (r_vector .* sin(phi_rad));
F_tip  = (2/pi()) * acos(exp(-f_tip));
% Secondary variables: Hub loss factor
f_hub  = (B/2   ) * (r_vector-R_root) ./ (r_vector .* sin(phi_rad));
F_hub  = (2/pi()) * acos(exp(-f_hub));
% Secondary variables: Global loss factor
if bool_loss == true
    F = F_tip .* F_hub;
else
    F = ones(size(mu_vector));
end

% Secondary variables: Actuator annuli thrust coefficient
if bool_glauert == true
    % Make branches and threshold
    Ct_axi_buhl_segment1    = 4 * a_axi .* (1-a_axi);                                                               % For low   a (such that CT <= 0/96F)
    Ct_axi_buhl_segment2    = 8/9 + (4*F - 40/9).* a_axi + (50/9 - 4*F) .* a_axi.^2;                                % For large a (such that CT >  0/96F)
    a_axi_threshold         = ( 18*F - 20 - 3*sqrt( (0.96*F) .* (50 - 36*F) + (12*F) .* (3*F - 4)) ) ./ (36*F - 50);% Is always 0.4
    % Now combine into result
    C_axi_rotor_annuli                          = Ct_axi_buhl_segment1;                                             % First copy classic expression into all stances
    C_axi_rotor_annuli(a_axi > a_axi_threshold) = Ct_axi_buhl_segment2(a_axi > a_axi_threshold);                    % Then replace stances where Ct exceeds 0.96 (which is equivalent to a>0.4)
else
    C_axi_rotor_annuli      = 4 * a_axi .* (1-a_axi);
end

% Secondary variables: Actuator annuli thrust
F_axi_rotor_annuli   = 0.5 * rho_inf * U_inf.^2 .* bin_area_vector .* C_axi_rotor_annuli .* F;

% Secondary variables: Actuator annuli moment
Q_tan_rotor_annuli   = 0.5 * rho_inf * U_inf.^2 .* (4 * mu_vector * lambda) .* (a_tan .* (1 - a_axi)) .* r_vector .* bin_area_vector .* F;

% Res
RES_axi = F_axi_rotor_annuli - F_axi_rotor_segments;
RES_tan = Q_tan_rotor_annuli - Q_tan_rotor_segments;

% Update
if bool_axi_solve == true
    a_axi = a_axi - 0.00001 * RES_axi;
end
if bool_tan_solve == true
    a_tan = a_tan - 0.00001 * RES_tan;
end

end

P_rotor_annuli   = sum(Q_tan_rotor_annuli  ) * Omega;
P_rotor_segments = sum(Q_tan_rotor_segments) * Omega;

CP_rotor_annuli    = P_rotor_annuli   / (0.5 * rho_inf * (pi*R^2) * U_inf^3);
CP_rotor_segments  = P_rotor_segments / (0.5 * rho_inf * (pi*R^2) * U_inf^3);
alpha_deg

[RES_axi ; RES_tan]
lambda
CP_rotor_annuli
CP_rotor_segments


