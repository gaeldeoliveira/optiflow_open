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
BC.R          = RWT.R     ;     % [m     ]  Rotor tip  radius (corresponds to R)
BC.R_root     = RWT.R_hub ;     % [m     ]  Rotor root radius
BC.Omega      = Omega_ref ;     % [rad/s ]  Rotational spped
BC.U_inf      = U_ref     ;     % [m/s   ]  Inflow speed
BC.B          = 3         ;     % [#     ]  Number of blades
BC.nu_inf     = 1.46e-5   ;     % [m2/s  ]  Fluid Viscosity (nu = mu/rho)
BC.rho_inf    = 1.225     ;     % [kg/m3 ]  Fluid Density
BC.N_bins     = 60        ;     % [#     ]  Number of computational elements along blade
% Parse inputs: localized induction / grid specifications
BC.N_az_bins  = 2+1      ;
% Parse inputs: inhomogeneous inflow specifications
BC.z0         = 0.0       ;    % [m     ]  Roughness lenght
BC.z_hub      = 120       ;    % [m     ]  Hub height (at which U_inf is defined)
BC.U_ref_fun                  = @(z  ) U_inf .* log(z/z0) ./ log(z_hub/z0);     % [m/s  ] Velocity at height z, as given by (nat) log law defined from roughness lenght, hub height, and unperturbed at speed at hub height
% Functional Inputs
%BC.c_over_R_of_mu_fun        = @(mu ) 0.1 * ones(size(mu));                                       % Make a constant   chord distribution
BC.c_over_R_of_mu_fun         = @(mu ) interp1(RWT.r_chord/ RWT.R, RWT.val_chord / RWT.R, mu);      % Chord distribution
%BC.theta_deg_of_mu_fun       = @(mu ) atan2(BC.U_inf, 70 * (mu*BC.R)) * 180 /pi();% - 5 ;         % Make a helical    twist distribution (don't reinterpolate)
BC.theta_deg_of_mu_fun        = @(mu ) interp1(RWT.r_twist_deg/ RWT.R, RWT.val_twist_deg      , mu);      % Twist distribution
BC.tc_of_mu_fun               = @(mu ) interp1(RWT.r_thickness/ RWT.R, RWT.val_thickness / 100, mu);      % Thickness distribution
BC.cl_of_alpha_deg_and_re_fun = @(alpha, Re, tc) polar_tensor_interpolator(polar_tensors, polar_tensors.cl_tensor, alpha, Re, tc);
BC.cd_of_alpha_deg_and_re_fun = @(alpha, Re, tc) polar_tensor_interpolator(polar_tensors, polar_tensors.cd_tensor, alpha, Re, tc);

% Solution details
BC.bool_loss       = true ;   % [bool.]  If true, activates tip loss correction (vanilla Prandtl)
BC.bool_glauert    = true  ;  % [bool.]  If true, activates Glauert  correction (vanilla Buhl   )
BC.bool_tan_solve  = true  ;  % [bool.]  If true, activates azimuthal induction (should be true )

BC.bool_prob_coef  = false ;
BC.bool_exceedance = false ;


[RES_axi, RES_tan, BD, BR] = bogo2_alfa8_BEMsolver(BC);

disp(['BR.P_rotor_annuli    = ' num2str(BR.P_rotor_annuli   )]);
disp(['BR.P_rotor_segments  = ' num2str(BR.P_rotor_segments(:)')]);
disp(['BR.CP_rotor_annuli   = ' num2str(BR.CP_rotor_annuli  )]);
disp(['BR.CP_rotor_segments = ' num2str(BR.CP_rotor_segments(:)')]);
disp(['BR.CT_rotor_annuli   = ' num2str(BR.CT_rotor_annuli  )]);
disp(['BR.CT_rotor_segments = ' num2str(BR.CT_rotor_segments(:)')]);
disp(['BR.CQ_axi_RBM        = ' num2str(BR.CQ_axi_root_bending_moment(:)')]);
disp(['BR.CQ_tan_RBM        = ' num2str(BR.CQ_tan_root_bending_moment(:)')]);
disp(['BR.RES_axi           = ' num2str(max(max(BR.RES_axi)))]);
disp(['BR.RES_tan           = ' num2str(max(max(BR.RES_tan)))]);

% % Now make 
% [RES_axi ; RES_tan]
% BR.lambda
% CP_rotor_annuli
% CP_rotor_segments
toc

% % And now plot stuff
% Support function, surface plots
seam_array = @(open_polar_mesh_array) [open_polar_mesh_array; open_polar_mesh_array(1,:)];
figure(10)
surf(seam_array(BD.y),seam_array(BD.z),seam_array(BD.U_inf_distribution )); view(2); xlabel('x'); ylabel('z'); colormap; shading interp; colorbar; title('Incomming wind speed (m/s)')
figure(1)
surf(seam_array(BD.y),seam_array(BD.z),seam_array(BR.a_axi              )); view(2); xlabel('x'); ylabel('z'); colormap; shading interp; colorbar; title('Axial induction')
figure(2)
surf(seam_array(BD.y),seam_array(BD.z),seam_array(BR.U_axi              )); view(2); xlabel('x'); ylabel('z'); colormap; shading interp; colorbar; title('Axial speed over rotor (m/s)')
figure(3)
surf(seam_array(BD.y),seam_array(BD.z),seam_array(BR.phi_deg            )); view(2); xlabel('x'); ylabel('z'); colormap; shading interp; colorbar; title('Inflow angle over rotor (deg)')


% Curve plots: axial induction 
figure(101)
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.a_axi, 0.8, BD.eta_vector(:,1))); 
xlabel('\eta - Azimuth angle [deg]'); ylabel('a_{axi}'); grid on; title('Axial induction factor'); hold on;
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.a_axi, 0.6, BD.eta_vector(:,1))); 
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.a_axi, 0.4, BD.eta_vector(:,1))); 
legend('\mu=r/R = 0.8', '\mu=r/R = 0.6', '\mu=r/R = 0.4');
% Curve plots: axial flow velocity on rotor plane
figure(102)
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.U_axi, 0.8, BD.eta_vector(:,1))); 
xlabel('\eta - Azimuth angle [deg]'); ylabel('U_{axi} [m/s]'); grid on; title('Axial flow speed'); hold on;
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.U_axi, 0.6, BD.eta_vector(:,1))); 
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.U_axi, 0.4, BD.eta_vector(:,1))); 
legend('\mu=r/R = 0.8', '\mu=r/R = 0.6', '\mu=r/R = 0.4');
% Curve plots: axial flow velocity on rotor plane
figure(103)
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.phi_deg, 0.8, BD.eta_vector(:,1))); 
xlabel('\eta - Azimuth angle [deg]'); ylabel('\phi [deg]'); grid on; title('Inflow angle on rotor'); hold on;
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.phi_deg, 0.6, BD.eta_vector(:,1))); 
plot(BD.eta_vector(:,1)*180/pi(), interp2(BD.mu_vector, BD.eta_vector, BR.phi_deg, 0.4, BD.eta_vector(:,1))); 
legend('\mu=r/R = 0.8', '\mu=r/R = 0.6', '\mu=r/R = 0.4');


%% Solver driver function
function [RES_axi, RES_tan, BD, BR] = bogo2_alfa8_BEMsolver(BC)
        % % Make discretization
        [~, ~, BD, ~] = bogo2_alfa8_RES_fun(BC, [], [], []);
        % Extract initial guess
        a_serial_0 = serialize_axi_tan(BD.a_axi_0, BD.a_tan_0);
        % Make serial function wrapper
        res_serial_fun = @(a_serial) bogo2_alfa8_RES_serialfun(BC, BD, a_serial);
        
        % % Now solve with non-linear solution algorithm
        % Make options for solution algorithm
        options = optimoptions('fsolve');
        options.StepTolerance           = 1e-18;
        options.MaxFunctionEvaluations  = 1e5;
        options.FunctionTolerance       = 1e-9;
        options.OptimalityTolerance     = 1e-9;
        options.MaxIterations           = 1e5;
        options.Display                 = 'off';
        % Solve!
        a_serial                        = fsolve(res_serial_fun, a_serial_0, options);
        
        % % Documentation round
        % Part primary variables
        [a_axi, a_tan] = part_axi_tan(a_serial);
        % Reshape them
        a_axi = reshape(a_axi, BC.N_az_bins-1, BC.N_bins-1);
        a_tan = reshape(a_tan, BC.N_az_bins-1, BC.N_bins-1);
        % Run documentation round
        [RES_axi, RES_tan, BD, BR] = bogo2_alfa8_RES_fun(BC, BD, a_axi, a_tan);
end

%% Residual evaluator serialization function
function [RES_serial] = bogo2_alfa8_RES_serialfun(BC, BD, a_serial)
% Serializes bogo2_alfa6_RES_fun inputs and outputs
    
    % Part primary variables
    [a_axi, a_tan] = part_axi_tan(a_serial);
    
    % Reshape them
    a_axi = reshape(a_axi, BC.N_az_bins-1, BC.N_bins-1);
    a_tan = reshape(a_tan, BC.N_az_bins-1, BC.N_bins-1);
    
    % Run function
    [RES_axi, RES_tan] = bogo2_alfa8_RES_fun(BC, BD, a_axi, a_tan);
    
    % Serialize outputs
    RES_serial = serialize_axi_tan(RES_axi, RES_tan);

end

%% Residual evaluator serialization function
function [RES_axi, RES_tan, BD, BR] = bogo2_alfa8_RES_fun(BC, BD, a_axi, a_tan)

    % Parse inputs: general
    R         = BC.R;                     % [m     ]  Rotor tip  radius (corresponds to R)
    R_root    = BC.R_root;                % [m     ]  Rotor root radius
    Omega     = BC.Omega;                 % [rad/s ]  Rotational spped
    U_inf     = BC.U_inf;                 % [m/s   ]  Inflow speed (reference, at hub height)
    B         = BC.B;                     % [#     ]  Number of blades
    nu_inf    = BC.nu_inf;                % [m2/s  ]  Fluid Viscosity (nu = mu/rho)
    rho_inf   = BC.rho_inf;               % [kg/m3 ]  Fluid Density
    N_bins    = BC.N_bins;                % [#     ]  Number of discretization bins
    % Parse inputs: new stuff for localized induction / grid specifications
    N_az_bins = BC.N_az_bins;
    % Parse inputs: new stuff for inhomogeneous inflow specifications
    z0        = BC.z0;                    % [m     ]  Roughness lenght
    z_hub     = BC.z_hub;                 % [m     ]  Hub height (at which U_inf is defined)
    %U_ref_fun BC.U_ref_fun;              % [m/s  ] Velocity at height z, as given by (nat) log law defined from roughness lenght, hub height, and unperturbed at speed at hub height

    % Parse inputs: solution procedure
    bool_loss      = BC.bool_loss;      % [bool  ]
    bool_glauert   = BC.bool_glauert;   % [bool  ]
    bool_tan_solve = BC.bool_tan_solve; %#ok<NASGU> % [bool  ]
    bool_prob_coef = BC.bool_prob_coef; % [bool  ]
    if isfield(BC, 'bool_exceedance')
        bool_exceedance= BC.bool_exceedance;% [bool  ]
    else
        bool_exceedance= false;
    end
    if isfield(BC, 'n_sigma_exceedance')
        n_sigma_exceedance = BC.n_sigma_exceedance;%[#]
    end

    % Parse inputs: functional
    c_over_R_of_mu_fun              = BC.c_over_R_of_mu_fun;              % Function of mu=r/R defining chord (c/R) distribution
    theta_deg_of_mu_fun             = BC.theta_deg_of_mu_fun;             % Function of mu=r/R defining twist [deg] distribution
    tc_of_mu_fun                    = BC.tc_of_mu_fun;                    % Function of mu=r/R defining thickness (t/c) distribution
    cl_of_alpha_deg_and_re_fun      = BC.cl_of_alpha_deg_and_re_fun;      % Function of alpha [deg] and Re [adim] defining lift polars
    cd_of_alpha_deg_and_re_fun      = BC.cd_of_alpha_deg_and_re_fun;      % Function of alpha [deg] and Re [adim] defining drag polars
    if or(bool_exceedance,bool_prob_coef)
        cl_of_alpha_deg_and_re_prob_fun = BC.cl_of_alpha_deg_and_re_prob_fun; % Function of alpha [deg], Re [adim], tc [adim] and mu [adim] defining expected value lift polars (probabilistic)
        cd_of_alpha_deg_and_re_prob_fun = BC.cd_of_alpha_deg_and_re_prob_fun; % Function of alpha [deg], Re [adim], tc [adim] and mu [adim] defining expected value drag polars (probabilistic)
    end
    
    % Make discretization if it is undefined
    if isempty(BD)
        % Radial     discretization (original: axisymmetric annuli)
        % BD.mu_bin_edges_vector  = linspace(R_root/R,1,N_bins);                                                              % [adim ] Edges of each discretization element
        % BD.mu_vector            = 0.5 * (BD.mu_bin_edges_vector(1:(end-1)) + BD.mu_bin_edges_vector(2:end));                % [adim ] Centers of each discretization element
        % BD.bin_width_vector     = R    * (BD.mu_bin_edges_vector(2:end) - BD.mu_bin_edges_vector(1:(end-1)));               % [m    ] Width of each discretization element (blade segment)
        % BD.bin_area_vector      = pi() * R.^2 * (BD.mu_bin_edges_vector(2:end).^2 - BD.mu_bin_edges_vector(1:(end-1)).^2);  % [m2   ] Area  of each discretization element (rotor annuli )
        % Polar grid discretization (new     : local induction    )
        % mu  = r/R [adim.] is the radial    coordinate (0 is hub center, 1 is tip. No deflection at this stage)
        % eta = r/R [rad  ] is the azimuthal coordinate (0 is tip up, pointing to the sky, at noon! same as 2pi) 
        % Define edge vectors and spread them over grid corner nodes
        [BD.mu_bin_edges_vector, BD.eta_bin_edges_vector]  = meshgrid(linspace(R_root/R,1,N_bins), linspace(0,2*pi(),N_az_bins));                         % [adim ] Edges of each discretization element (i direction azimuthal, j direction 
        % Allocate arrays for cell centered data
        BD.mu_vector            = zeros(size(BD.mu_bin_edges_vector) - [1 1]);
        BD.eta_vector           = zeros(size(BD.mu_bin_edges_vector) - [1 1]);
        BD.bin_width_vector     = zeros(size(BD.mu_bin_edges_vector) - [1 1]);
        BD.bin_angle_vector     = zeros(size(BD.mu_bin_edges_vector) - [1 1]);
        BD.bin_area_vector      = zeros(size(BD.mu_bin_edges_vector) - [1 1]);
        % Fill static cell centered data from grid corner data
        for n_azi_bin = 1:(size(BD.mu_bin_edges_vector,1)-1)
            for n_rad_bin = 1:(size(BD.mu_bin_edges_vector,2)-1)
                % Radial    center of each grid cell                             [adim ] (r/R)
                BD.mu_vector(       n_azi_bin, n_rad_bin) = 0.5 * (BD.mu_bin_edges_vector(n_azi_bin   , n_rad_bin+1)    + BD.mu_bin_edges_vector( n_azi_bin, n_rad_bin)   );
                % Azimuthal center of each grid cell                             [rad  ] 
                BD.eta_vector(      n_azi_bin, n_rad_bin) = 0.5 * (BD.eta_bin_edges_vector(n_azi_bin+1, n_rad_bin  )    + BD.eta_bin_edges_vector(n_azi_bin, n_rad_bin)   );
                % Radial    width of each grid cell              (blade segment) [m    ] 
                BD.bin_width_vector(n_azi_bin, n_rad_bin) = R   * (BD.mu_bin_edges_vector(n_azi_bin   , n_rad_bin+1)    - BD.mu_bin_edges_vector( n_azi_bin, n_rad_bin)   );
                % Azimuthal width of each grid cell   (fraction of rotor annuli) [rad  ] 
                BD.bin_angle_vector(n_azi_bin, n_rad_bin) =       (BD.eta_bin_edges_vector(n_azi_bin+1, n_rad_bin  )    - BD.eta_bin_edges_vector(n_azi_bin, n_rad_bin)   );
                % Area  of each grid cell             (fraction of rotor annuli) [m2   ] 
                BD.bin_area_vector( n_azi_bin, n_rad_bin) = 0.5 * BD.bin_angle_vector(n_azi_bin, n_rad_bin) .* R.^2 .* ...
                                                                  (BD.mu_bin_edges_vector(n_azi_bin   , n_rad_bin+1).^2 - BD.mu_bin_edges_vector( n_azi_bin, n_rad_bin).^2);
            end
        end
        
        % Immatuable scalars (new! for inhomogeneous incoming flow velocity)
        BD.z                    = R * BD.mu_vector .* cos(BD.eta_vector) + z_hub;       % [m    ] Height from ground as function of radial (mu=r/R) and azimutal coordinate (eta in rad, 0 tip points to the sky)
        BD.y                    = R * BD.mu_vector .* sin(BD.eta_vector) + 0    ;       % [m    ] Side distance from hub center as function of radial (mu=r/R) and azimutal coordinate (eta in rad, 0 tip points to the sky)
        if z0 > 0    
            BD.U_inf_distribution   = U_inf .* log(BD.z/z0) ./ log(z_hub/z0);           % [m/s  ] Velocity at height z, as given by (nat) log law defined from roughness lenght, hub height, and unperturbed at speed at hub height
        else
            BD.U_inf_distribution   = U_inf .* ones(size(BD.z));                        % [m/s  ] Velocity at height z, as given by (nat) log law defined from roughness lenght, hub height, and unperturbed at speed at hub height
        end
        % Useful verification plot
        % surf(BD.y,BD.z,BD.U_inf_distribution); view(2); xlabel('x'); ylabel('z'); colormap; shading interp

        % Immutable scalars
        BD.lambda = Omega * R / U_inf;                                                  % [adim ] Tip speed ratio
        % Immutable arrays
        BD.r_vector             = BD.mu_vector * R;                                     % [m    ] Radius                     (at bin/cell center)
        BD.c_vector             = c_over_R_of_mu_fun(BD.mu_vector) * R;                 % [m    ] Chord                      (at bin/cell center)
        BD.tc_vector            = tc_of_mu_fun(BD.mu_vector);                           % [adim ] Thickness over chord ratio (at bin/cell center)
        BD.theta_deg            = theta_deg_of_mu_fun(BD.mu_vector);                    % [deg  ] Thickness over chord ratio (at bin/cell center)
        % Primary Variables (first guess/input)
        BD.a_axi_0 = zeros(size(BD.mu_vector));                                         % [adim ] Initial guess for axial      induction factor
        BD.a_tan_0 = zeros(size(BD.mu_vector));                                         % [adim ] Initial guess for tangential induction factor
        % Stop here
        BR      = [];
        RES_axi = [];
        RES_tan = [];
        return
    % Otherwise extract data and compute residuals
    else
        % Extract discretization
        mu_bin_edges_vector  = BD.mu_bin_edges_vector;  %#ok<NASGU> % [adim ] Edges of each grid cell (radial    information, mu=r/R , 0 means hub center)  
        eta_bin_edges_vector = BD.eta_bin_edges_vector; %#ok<NASGU> % [rad  ] Edges of each grid cell (azimuthal information, eta in radians, 0 points up)
        mu_vector            = BD.mu_vector;             % [adim ] Centers of each discretization element
        eta_vector           = BD.eta_vector;            %#ok<NASGU> % [rad  ] Azimuthal center of each grid cell
        bin_width_vector     = BD.bin_width_vector;      % [m    ] Width of each discretization element (blade segment)
        bin_angle_vector     = BD.bin_angle_vector;      % [rad  ] Azimuthal width of each grid cell   (fraction of rotor annuli) 
        bin_area_vector      = BD.bin_area_vector;       % [m2   ] Area  of each discretization element (fraction of rotor annuli )
        % Immatuable scalars (new! for inhomogeneous incoming flow velocity)
        z                    = BD.z;                     %#ok<NASGU> % [m    ] Height from ground as function of radial (mu=r/R) and azimutal coordinate (eta in rad, 0 tip points to the sky)
        y                    = BD.y;                     %#ok<NASGU> % [m    ] Side distance from hub center as function of radial (mu=r/R) and azimutal coordinate (eta in rad, 0 tip points to the sky)
        U_inf_distribution   = BD.U_inf_distribution;    % [m/s  ] Velocity at height z, as given by (nat) log law defined from roughness lenght, hub height, and unperturbed at speed at hub height
        
        % Extract immutable scalars
        lambda               = BD.lambda;                % [adim ] Tip speed ratio
        % Extract immutable arrays 
        r_vector             = BD.r_vector;
        c_vector             = BD.c_vector;
        tc_vector            = BD.tc_vector;
        theta_deg            = BD.theta_deg;
        
        % Secondary variables: Speeds
        U_axi = U_inf_distribution .* (1 - a_axi);
        U_tan = Omega * r_vector   .* (1 + a_tan);
        W     = sqrt(U_axi.^2 + U_tan.^2);
        
        % Secondary variables: Reynolds number
        Re_vector = W .* c_vector / nu_inf;
        
        % Secondary variables: Angles
        phi_rad   = atan2(U_axi, U_tan);
        phi_deg   = phi_rad * 180 /pi();
        
        alpha_deg = phi_deg - theta_deg;
        
        % Secondary variables: Airfoil coefficients
        if bool_prob_coef == false
            cl_vector            = cl_of_alpha_deg_and_re_fun(alpha_deg, Re_vector, tc_vector);
            cd_vector            = cd_of_alpha_deg_and_re_fun(alpha_deg, Re_vector, tc_vector);
        else
            cl_vector            = cl_of_alpha_deg_and_re_prob_fun(alpha_deg, Re_vector, tc_vector, mu_vector);
            cd_vector            = cd_of_alpha_deg_and_re_prob_fun(alpha_deg, Re_vector, tc_vector, mu_vector);
            if bool_exceedance == true
                cl_vector        = cl_of_alpha_deg_and_re_fun(alpha_deg + n_sigma_exceedance * mu_vector, Re_vector, tc_vector);
                cd_vector        = cl_of_alpha_deg_and_re_fun(alpha_deg + n_sigma_exceedance * mu_vector, Re_vector, tc_vector);
            end
        end
        
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
            Ct_axi_buhl_segment1    = 4 * a_axi .* F .* (1-a_axi);                                                           % For low   a (such that CT <= 0/96F) (this way to apply F (4*a*(1-a)*F) has been checked: it is consistent with the upper streak (despite being different from the expression of Gijs' Chap. 8, which applies to a different context)
            Ct_axi_buhl_segment2    = 8/9 + (4*F - 40/9).* a_axi + (50/9 - 4*F) .* a_axi.^2;                                 % For large a (such that CT >  0/96F) (this is the vanilla expression of Buhl, is know to work well)
            a_axi_threshold         = ( 18*F - 20 - 3*sqrt( (0.96*F) .* (50 - 36*F) + (12*F) .* (3*F - 4)) ) ./ (36*F - 50); % Is always 0.4
            % Now combine into result
            C_axi_rotor_annuli                          = Ct_axi_buhl_segment1;                                              % First copy classic expression into all stances
            C_axi_rotor_annuli(a_axi > a_axi_threshold) = Ct_axi_buhl_segment2(a_axi > a_axi_threshold);                     % Then replace stances where Ct exceeds 0.96 (which is equivalent to a>0.4)
        else
            C_axi_rotor_annuli      = 4 * a_axi .* F .* (1-a_axi);
        end
        
        % Secondary variables: Actuator annuli thrust
        F_axi_rotor_annuli   = 0.5 * rho_inf .* U_inf_distribution.^2 .* bin_area_vector .* C_axi_rotor_annuli;                                                % Removed F from here, because it is now everywhere above. Version used for PhD computations had F here and not above, except for high CT where it was also above. That was formally incorrect but had strictly no effect on results because all points of interest must have had CT<0.96.
        
        % Secondary variables: Actuator annuli moment
        Q_tan_rotor_annuli   = 0.5 * rho_inf .* U_inf_distribution.^2 .* (4 * mu_vector * lambda) .* (a_tan .* (1 - a_axi)) .* r_vector .* bin_area_vector;    % TO CHECK: Removed F from here. Phd Computation had it. Difference on Cp of innwind blade: .48798 (without F) instead of 0.48814 (with F, as in PhD computations). 0.0314% power coefficient suggest this should be checked, but with modest priority.
        
        % Residuals (consistency now (localized induction) means equalization of pressure per surface area, hence scalling of F_axi_rotor_annuli and Q_tan_rotor_annuli)
        RES_axi = F_axi_rotor_annuli .* (2*pi()./bin_angle_vector) - F_axi_rotor_segments;
        RES_tan = Q_tan_rotor_annuli .* (2*pi()./bin_angle_vector) - Q_tan_rotor_segments;
    end
    
    
    % % If this is the final round
    if nargout > 2
        % Postprocess: Power
        P_rotor_annuli       = sum(sum(Q_tan_rotor_annuli))    * Omega;
        P_rotor_segments     =     sum(Q_tan_rotor_segments,2) * Omega;
        % Postprocess: Power Coefficient
        CP_rotor_annuli      = P_rotor_annuli   / (0.5 * rho_inf * (pi*R^2) * U_inf^3);                           %#ok<NASGU>
        CP_rotor_segments    = P_rotor_segments / (0.5 * rho_inf * (pi*R^2) * U_inf^3);                           %#ok<NASGU>
        % Postprocess: Thrust
        T_rotor_annuli       = sum(sum(F_axi_rotor_annuli));
        T_rotor_segments     =     sum(F_axi_rotor_segments,2);
        % Postprocess: Thrust Coefficient
        CT_rotor_annuli      = T_rotor_annuli   / (0.5 * rho_inf * (pi*R^2) * U_inf^2);                           %#ok<NASGU>
        CT_rotor_segments    = T_rotor_segments / (0.5 * rho_inf * (pi*R^2) * U_inf^2);                           %#ok<NASGU>
        % Postprocess: Out of plane bending moment
        Q_axi_rotor_segments = F_axi_rotor_segments .* r_vector;
        % Postprocess: Root bending moments
        Q_axi_root_bending_moment  = sum(Q_axi_rotor_segments,2);
        Q_tan_root_bending_moment  = sum(Q_tan_rotor_segments,2);
        % Postprocess: Root bending moment coefficients
        CQ_axi_root_bending_moment = Q_axi_root_bending_moment / (0.5 * rho_inf * (pi*R^2) * U_inf^2 * R);     %#ok<NASGU>
        CQ_tan_root_bending_moment = Q_tan_root_bending_moment / (0.5 * rho_inf * (pi*R^2) * U_inf^2 * R);     %#ok<NASGU>

        % % And output all variables
        c = who();
        BR = struct();
        for n_var = 1:length(c)
            BR.(c{n_var}) = eval(c{n_var});
        end
    end
end

