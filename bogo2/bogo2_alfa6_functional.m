% A very simple BEM code bogo2_alfa6
addpath src

% % Define BemCase
% Scalar Inputs
BC.R       = 1         ;     % [m     ]  Rotor tip  radius (corresponds to R)
BC.R_root  = 0.2       ;     % [m     ]  Rotor root radius
BC.Omega   = 70      ;     % [rad/s ]  Rotational spped
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
BC.theta_deg_of_mu_fun        = @(mu ) atan2(BC.U_inf, 70 * (mu*BC.R)) * 180 /pi() - 5 ;                        % Make a helical    twist distribution (don't reinterpolate)
BC.cl_of_alpha_deg_and_re_fun = @(alpha_deg, Re) 2*pi*sin(alpha_deg*pi/180).*cos(alpha_deg*pi/180);
BC.cd_of_alpha_deg_and_re_fun = @(alpha_deg, Re) 0.01 * (1+cl_of_alpha_deg_and_re_fun(alpha_deg, Re)).^2;
BC.cd_of_alpha_deg_and_re_fun = @(alpha_deg, Re) zeros(size(alpha_deg));


% % Make discretization
[~, ~, BD, ~] = bogo2_alfa6_RES_fun(BC, [], [], []);

% % Extract initial guess
a_axi = BD.a_axi_0;
a_tan = BD.a_tan_0;

% % Solve with simple iterative loop
for n = 1:100000
    % Secondary variables: Positions and Sizes
    [RES_axi, RES_tan] = bogo2_alfa6_RES_fun(BC, BD, a_axi, a_tan);
    % Update
    if BC.bool_axi_solve == true
        % a_axi = a_axi - 0.0001 * RES_axi;
        a_axi = a_axi - min(0.0001 * RES_axi, 0.001*ones(size(RES_axi)));
    end
    if BC.bool_tan_solve == true
        % a_tan = a_tan - 0.0001 * RES_tan;
        a_tan = a_tan - min(0.0001 * RES_tan, 0.001*ones(size(RES_tan)));
    end 
end
[RES_axi, RES_tan, BD, BR] = bogo2_alfa6_RES_fun(BC, BD, a_axi, a_tan);

% % % Now solve with simple iterative with serialized interface
% % Serialize initial guess
a_serial_0 = serialize_axi_tan(BD.a_axi_0, BD.a_tan_0);
% % Make serial function
res_serial_fun = @(a_serial) bogo2_alfa6_RES_serialfun(BC, BD, a_serial);
% % Start up
% a_serial = a_serial_0;
% % Iterate
% for n = 1:10000
%     a_serial = a_serial - 0.001 * res_serial_fun(a_serial);
% end

% % Now solve with non-linear solution algorithm
% Serialize initial guess
a_serial_0 = serialize_axi_tan(BD.a_axi_0, BD.a_tan_0);
% Make serial function
res_serial_fun = @(a_serial) bogo2_alfa6_RES_serialfun(BC, BD, a_serial);
% Make options for solution algorithm
options = optimoptions('fsolve');
options.StepTolerance           = 1e-18;
options.MaxFunctionEvaluations  = 1e5;
options.FunctionTolerance       = 1e-9;
options.OptimalityTolerance     = 1e-9;
options.MaxIterations = 1e5
% % Solve!
a_serial       = fsolve(res_serial_fun, a_serial_0, options);

% % Make last round for documentation
[a_axi2, a_tan2] = part_axi_tan(a_serial);
[RES_axi2, RES_tan2, BD2, BR2] = bogo2_alfa6_RES_fun(BC, BD, a_axi2, a_tan2);


disp(['BR.CP_rotor_annuli   = ' num2str(BR.CP_rotor_annuli  )])
disp(['BR.CP_rotor_segments = ' num2str(BR.CP_rotor_segments)])
disp(['BR.RES_axi           = ' num2str(max(BR.RES_axi)     )])
disp(['BR.RES_tan           = ' num2str(max(BR.RES_tan)     )])

disp(['BR2.CP_rotor_annuli   = ' num2str(BR2.CP_rotor_annuli  )])
disp(['BR2.CP_rotor_segments = ' num2str(BR2.CP_rotor_segments)])
disp(['BR2.RES_axi           = ' num2str(max(BR2.RES_axi)     )])
disp(['BR2.RES_tan           = ' num2str(max(BR2.RES_tan)     )])

% % Now make 
% [RES_axi ; RES_tan]
% BR.lambda
% CP_rotor_annuli
% CP_rotor_segments


