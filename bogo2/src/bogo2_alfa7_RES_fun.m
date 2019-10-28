function [RES_axi, RES_tan, BD, BR] = bogo2_alfa7_RES_fun(BC, BD, a_axi, a_tan)

    % Parse inputs
    R       = BC.R;                     % [m     ]  Rotor tip  radius (corresponds to R)
    R_root  = BC.R_root;                % [m     ]  Rotor root radius
    Omega   = BC.Omega;                 % [rad/s ]  Rotational spped
    U_inf   = BC.U_inf;                 % [m/s   ]  Inflow speed
    B       = BC.B;                     % [#     ]  Number of blades
    nu_inf  = BC.nu_inf;                % [m2/s  ]  Fluid Viscosity (nu = mu/rho)
    rho_inf = BC.rho_inf;               % [kg/m3 ]  Fluid Density
    N_bins  = BC.N_bins;                % [#     ]  Number of discretization bins

    % Solution procedure
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

    % Functional Inputs
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
        % Discretization
        BD.mu_bin_edges_vector  = linspace(R_root/R,1,N_bins);                                                              % [adim ] Edges of each discretization element
        BD.mu_vector            = 0.5 * (BD.mu_bin_edges_vector(1:(end-1)) + BD.mu_bin_edges_vector(2:end));                % [adim ] Centers of each discretization element
        BD.bin_width_vector     = R    * (BD.mu_bin_edges_vector(2:end) - BD.mu_bin_edges_vector(1:(end-1)));               % [m    ] Width of each discretization element (blade segment)
        BD.bin_area_vector      = pi() * R.^2 * (BD.mu_bin_edges_vector(2:end).^2 - BD.mu_bin_edges_vector(1:(end-1)).^2);  % [m2   ] Area  of each discretization element (rotor annuli )
        % Immutable scalars
        BD.lambda = Omega * R / U_inf;                                                                                      % [adim ] Tip speed ratio
        % Immutable arrays
        BD.r_vector             = BD.mu_vector * R;
        BD.c_vector             = c_over_R_of_mu_fun(BD.mu_vector) * R;
        BD.tc_vector            = tc_of_mu_fun(BD.mu_vector);
        BD.theta_deg            = theta_deg_of_mu_fun(BD.mu_vector);
        % Primary Variables (first guess/input)
        BD.a_axi_0 = zeros(size(BD.mu_vector));                                                                             % [adim ] Initial guess for axial      induction factor
        BD.a_tan_0 = zeros(size(BD.mu_vector));                                                                             % [adim ] Initial guess for tangential induction factor
        % Stop here
        BR      = [];
        RES_axi = [];
        RES_tan = [];
        return
    % Otherwise extract data and compute residuals
    else
        % Extract discretization
        mu_bin_edges_vector = BD.mu_bin_edges_vector;   % [adim ] Edges of each discretization element
        mu_vector           = BD.mu_vector;             % [adim ] Centers of each discretization element
        bin_width_vector    = BD.bin_width_vector;      % [m    ] Width of each discretization element (blade segment)
        bin_area_vector     = BD.bin_area_vector;       % [m2   ] Area  of each discretization element (rotor annuli )
        % Extract immutable scalars
        lambda              = BD.lambda;                % [adim ] Tip speed ratio
        % Extract immutable arrays 
        r_vector            = BD.r_vector;
        c_vector            = BD.c_vector;
        tc_vector           = BD.tc_vector;
        theta_deg           = BD.theta_deg;
        
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
            Ct_axi_buhl_segment1    = 4 * a_axi .*F .* (1-a_axi);                                                           % For low   a (such that CT <= 0/96F) (this way to apply F (4*a*(1-a)*F) has been checked: it is consistent with the upper streak (despite being different from the expression of Gijs' Chap. 8, which applies to a different context)
            Ct_axi_buhl_segment2    = 8/9 + (4*F - 40/9).* a_axi + (50/9 - 4*F) .* a_axi.^2;                                % For large a (such that CT >  0/96F)
            a_axi_threshold         = ( 18*F - 20 - 3*sqrt( (0.96*F) .* (50 - 36*F) + (12*F) .* (3*F - 4)) ) ./ (36*F - 50);% Is always 0.4
            % Now combine into result
            C_axi_rotor_annuli                          = Ct_axi_buhl_segment1;                                             % First copy classic expression into all stances
            C_axi_rotor_annuli(a_axi > a_axi_threshold) = Ct_axi_buhl_segment2(a_axi > a_axi_threshold);                    % Then replace stances where Ct exceeds 0.96 (which is equivalent to a>0.4)
        else
            C_axi_rotor_annuli      = 4 * a_axi .*F .* (1-a_axi);
        end
        
        % Secondary variables: Actuator annuli thrust
        F_axi_rotor_annuli   = 0.5 * rho_inf * U_inf.^2 .* bin_area_vector .* C_axi_rotor_annuli;                           % Removed F from here, because it is now everywhere above. Version used for PhD computations had F here and not above, except for high CT where it was also above. That was formally incorrect but had strictly no effect on results because all points of interest must have had CT<0.96.
        
        % Secondary variables: Actuator annuli moment
        Q_tan_rotor_annuli   = 0.5 * rho_inf * U_inf.^2 .* (4 * mu_vector * lambda) .* (a_tan .* (1 - a_axi)) .* r_vector .* bin_area_vector;      % TO CHECK: Removed F from here. Phd Computation had it. Difference on Cp of innwind blade: .48798 (without F) instead of 0.48814 (with F, as in PhD computations). 0.0314% power coefficient suggest this should be checked, but with modest priority.
        
        % Residuals
        RES_axi = F_axi_rotor_annuli - F_axi_rotor_segments;
        RES_tan = Q_tan_rotor_annuli - Q_tan_rotor_segments;
    end
    
    
    % % If this is the final round
    if nargout > 2
        % Postprocess: Power
        P_rotor_annuli   = sum(Q_tan_rotor_annuli  ) * Omega;
        P_rotor_segments = sum(Q_tan_rotor_segments) * Omega;
        % Postprocess: Power Coefficient
        CP_rotor_annuli    = P_rotor_annuli   / (0.5 * rho_inf * (pi*R^2) * U_inf^3);                           %#ok<NASGU>
        CP_rotor_segments  = P_rotor_segments / (0.5 * rho_inf * (pi*R^2) * U_inf^3);                           %#ok<NASGU>
        % Postprocess: Thrust
        T_rotor_segments   = sum(F_axi_rotor_segments);                               
        T_rotor_annuli     = sum(F_axi_rotor_annuli);                                 
        % Postprocess: Thrust Coefficient
        CT_rotor_annuli    = T_rotor_annuli   / (0.5 * rho_inf * (pi*R^2) * U_inf^2);                           %#ok<NASGU>
        CT_rotor_segments  = T_rotor_segments / (0.5 * rho_inf * (pi*R^2) * U_inf^2);                           %#ok<NASGU>
        % Postprocess: Out of plane bending moment
        Q_axi_rotor_segments = F_axi_rotor_segments .* r_vector;
        % Postprocess: Root bending moments
        Q_axi_root_bending_moment  = sum(Q_axi_rotor_segments);                       
        Q_tan_root_bending_moment  = sum(Q_tan_rotor_segments);                       
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
