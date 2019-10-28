% % Design problem settings (common to study)
PB.make_optimization        = true  ;
PB.MaxGenerations           =  60   ;
PB.PopulationSize           = 180   ;
PB.ParetoFraction           = 0.4   ;
PB.UseParallel              =  true ;
PB.base_blade_geometry_file = '../rotor_integration/graph_digitization/inwind_moded_RWT.mat'  ;
PB.std_distribution_file    = '../rotor_integration/probabilistic/std_distribution_scaled.mat';
if PB.make_optimization == true
    PB.start_datestr            = ['DO_' , datestr(now,30)];
else
    PB.start_datestr            = ['DT_' , datestr(now,30)];
end

% % Make folder to store results
mkdir([PB.start_datestr])

% % Stat, FFA foils
% Set
PB.polar_tensors_file       = '../rotor_integration/airfoil_families/FFA/FFA_freeC.mat'     ;
PB.bool_prob_coef           = false;
PB.case_id                  = 'stat_FFA_free_limit_case';

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

% % Modify polar tensors for limit case
% Define ideal cl fucntion
cl0          = 0.3867;
cl_ideal_fun = @(alpha) 2*pi*sin(alpha*pi/180).*cos(alpha*pi/180) + cl0;
% Make backup
polar_tensors_backup = polar_tensors;
% Make shorthand for range of angles of attack
alpha_range = polar_tensors.alpha_range;
% Modify range of polar_tensors
for n_alpha = 1:length(polar_tensors.alpha_range)
    for n_Re = 1:length(polar_tensors.Re_range)
        n_tc = 1; % for n_tc = 1 (%corresponds to lowest thickness, 21% mismatched to 24%)
            polar_tensors.cl_tensor(n_alpha, n_Re,n_tc) = cl_ideal_fun(polar_tensors.alpha_range(n_alpha));
            polar_tensors.cd_tensor(n_alpha, n_Re,n_tc) = 0;    % Ideal case is 0 drag
        % end
    end
end

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

% Solution details
BC.bool_loss            = true;
BC.bool_glauert         = true;
BC.bool_tan_solve       = true;
BC.bool_prob_coef       = PB.bool_prob_coef;
BC.bool_exceedance      = false; % for computations, switched on at end!

% % Anounce program entry
disp(['======================================================']); %#ok<NBRAK>
disp(['==== BOGO2 with GAM ==================================']); %#ok<NBRAK>
disp(['============= Welcome ================================']); %#ok<NBRAK>
disp(['============================== Gael de Oliveira ======']); %#ok<NBRAK>
disp(['======================================================']); %#ok<NBRAK>
disp(['Base geometry   file : ' , PB.base_blade_geometry_file]);
disp(['Polar tensors   file : ' , PB.polar_tensors_file]);
disp(['AOA probability file : ' , PB.std_distribution_file]);
disp(['==== Compute reference case : ' , datestr(now),  ' ===']);

% % Start with solution on reference blade
[RES_axi, RES_tan, BD, BR] = bogo2_alfa7_BEMsolver(BC);

disp(['==== Done with computation  : ' , datestr(now),  ' ===']);
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

if PB.make_optimization == true
    % % Now make optimization
    % [RES_axi ; RES_tan]
    % BR.lambda
    % CP_rotor_annuli
    % CP_rotor_segments
    
    
    % Make reference case design vector
    A_free_0 = ones(1,8);
    % Test cost function
    [cf] = modify_and_solve_BEM_case_serialfun(BC, A_free_0);
    
    % Make upper and lower bounds
    A_free_lb = ones(1,8) * 0.5;
    A_free_ub = ones(1,8) * 1.5;
    
    % Make genetic algorithm options structure
    options   = optimoptions('gamultiobj');
    options.MaxGenerations = PB.MaxGenerations;
    options.PopulationSize = PB.PopulationSize;
    options.ParetoFraction = PB.ParetoFraction;
    options.UseParallel    = PB.UseParallel   ;
    options.PlotFcn        = @gaplotpareto    ;
    options.Display        = 'iter'           ;
    % Start parallel pool of applicable
    if PB.UseParallel == true
        try
            parpool('local');
        catch
            disp('BOGO2 with GAM : Parallel pool unavailable');
        end
    end
    
    % Make function wrapper
    cf_wrapper_fun = @(A_free) modify_and_solve_BEM_case_serialfun(BC, A_free);
    
    % Test it on all relevant cases
    disp(['==== Prepare for Optimization : ' , datestr(now), ' ====']);
    disp(['Reference   case (CP, CQ_axi_RBM) : ' , num2str(cf_wrapper_fun(A_free_0))]);
    disp(['Upper bound case (CP, CQ_axi_RBM) : ' , num2str(cf_wrapper_fun(A_free_lb))]);
    disp(['Lower bound case (CP, CQ_axi_RBM) : ' , num2str(cf_wrapper_fun(A_free_ub))]);
    
    % Now announce optimization
    disp(['==== Launch GAM  Optimization : ' , datestr(now), ' =']);
    disp(['Max generations : ' num2str(options.MaxGenerations) ]);
    disp(['Population size : ' num2str(options.PopulationSize) ]);
    disp(['Pareto fraction : ' num2str(options.ParetoFraction) ]);
    % Set rng state
    rng('default')
    % And launch it!
    [A_free_opt,fval,exitFlag,output,population,scores] = gamultiobj(cf_wrapper_fun,8,[],[],[],[],A_free_lb,A_free_ub,[],options);
    % When done, announce completion
    disp(['==== Optimization reached end : ' , datestr(now), ' ====']);
    % Postrocess: first allocate space
    CP_cell         = cell(size(A_free_opt,1), 1);
    CQ_axi_RBM_cell = cell(size(A_free_opt,1), 1);
    A_chord_cell    = cell(size(A_free_opt,1), 1);
    A_theta_cell    = cell(size(A_free_opt,1), 1);
    BC_moded_cell   = cell(size(A_free_opt,1), 1);
    BD_cell         = cell(size(A_free_opt,1), 1);
    BR_cell         = cell(size(A_free_opt,1), 1);
    % Then fill it
    for n_opt = 1:size(A_free_opt,1)
        [CP_cell{n_opt}, CQ_axi_RBM_cell{n_opt}, A_chord_cell{n_opt}, A_theta_cell{n_opt}, BC_moded_cell{n_opt}, BD_cell{n_opt}, BR_cell{n_opt}] = modify_and_solve_BEM_case(BC, A_free_opt(n_opt,:));
    end
else
    disp(['==== No optimization requested =======================']);  %#ok<NBRAK>
end

% Save results
save([PB.start_datestr, '/', PB.case_id, '.mat']);
    
disp(['==== All variables were saved : ' , datestr(now), ' =']);
disp(['======================================================']); %#ok<NBRAK>
disp(['==== BOGO2 with GAM ==================================']); %#ok<NBRAK>
disp(['============= See you soon ! =========================']); %#ok<NBRAK>
disp(['============================== Gael de Oliveira ======']); %#ok<NBRAK>
disp(['======================================================']); %#ok<NBRAK>

if PB.make_optimization == true
    % % Post stuff
    % Plot pareto front
    figure;
    plot(cf(1), cf(2), 'o'); hold on; grid on; plot(fval(:,1), fval(:,2), 'x')
    % Plot effect aoa perturbations on polar curve of expected values
    alpha_range = -5:0.2:20;
    figure;
    subplot(2,2,1)
    plot(alpha_range, BC.cl_of_alpha_deg_and_re_fun(     alpha_range,9e6,0.21  )); hold on; xlabel('\alpha [deg]')
    plot(alpha_range, BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range,9e6,0.21,1)); grid on; ylabel('C_l')
    subplot(2,2,2)
    plot(alpha_range, BC.cd_of_alpha_deg_and_re_fun(     alpha_range,9e6,0.21  )); hold on; xlabel('\alpha [deg]')
    plot(alpha_range, BC.cd_of_alpha_deg_and_re_prob_fun(alpha_range,9e6,0.21,1)); grid on; ylabel('C_d')
    subplot(2,2,3)
    plot(alpha_range, BC.cl_of_alpha_deg_and_re_fun(     alpha_range,9e6,0.21  ) ./ BC.cd_of_alpha_deg_and_re_fun(     alpha_range,3e6,0.21  )); hold on; xlabel('\alpha [deg]')
    plot(alpha_range, BC.cl_of_alpha_deg_and_re_prob_fun(alpha_range,9e6,0.21,1) ./ BC.cd_of_alpha_deg_and_re_prob_fun(alpha_range,3e6,0.21,1)); grid on; ylabel('L/D')
end


function [cf] = modify_and_solve_BEM_case_serialfun(BC_base, A_free)
    % Make modified BEM case and solve it
    [CP, CQ_axi_RBM, A_chord, A_theta, BC_moded, BD, BR] = modify_and_solve_BEM_case(BC_base, A_free); %#ok<ASGLU>
    % Serialize outputs (and invert for minimization as cost function if needed)
    cf = zeros(1,2);
    cf(1) =   CP        ;
    cf(2) =   CQ_axi_RBM;
    % Eliminate case if convergence was not good enough
    if max(abs(BR.RES_axi)) + max(abs(BR.RES_tan)) > 1e-6
        cf(1) = 0;
        cf(2) = Inf;
    end
end

function [CP, CQ_axi_RBM, A_chord, A_theta, BC_moded, BD, BR] = modify_and_solve_BEM_case(BC_base, A_free)
    % Generate modify shapes form A_free (first half on chord distribution,
    % second half on twist distribution)
    A_chord = ones(1,8); A_chord(5:8) = A_free(1:4);
    A_theta = ones(1,8); A_theta(5:8) = A_free(5:8);
    
    % Make modified BEM case
    [BC_moded] = modify_BEMcase(BC_base, A_chord, A_theta);
    
%     % Old
%     % Solve it
%     [RES_axi, RES_tan, BD, BR] = bogo2_alfa7_BEMsolver(BC_moded); %#ok<ASGLU>
%         
%     % Extract outputs
%     CP          = BR.CP_rotor_segments         ;
%     CQ_axi_RBM  = BR.CQ_axi_root_bending_moment;

    % Solve it
    if BC_moded.bool_prob_coef == false
        % Once for CP and CQ in static mode
        [RES_axi, RES_tan, BD, BR] = bogo2_alfa7_BEMsolver(BC_moded); %#ok<ASGLU>
        
        % Extract outputs
        CP          = BR.CP_rotor_segments         ;
        CQ_axi_RBM  = BR.CQ_axi_root_bending_moment;
    else
        % Solve in prob mode for CP
        BC_moded.bool_exceedance = false;
        [RES_axi, RES_tan, BD, BR] = bogo2_alfa7_BEMsolver(BC_moded); %#ok<ASGLU>
        % Extract CP
        CP          = BR.CP_rotor_segments              ;
        % Post-process to get coherent gust load
        % Add gusts to angle of attack but don't change inflow
        BC_moded.bool_exceedance = true;
        % Set magnitude of gust and recompute forces
        BC_moded.n_sigma_exceedance = 0.0;
        [RES_axi_exc_00, RES_tan_exc_00, BD_exc_00, BR_exc_00] = bogo2_alfa7_RES_fun(BC_moded, BD, BR.a_axi, BR.a_tan); %#ok<ASGLU>
        % Set magnitude of gust and recompute forces
        BC_moded.n_sigma_exceedance = 0.5;
        [RES_axi_exc_05, RES_tan_exc_05, BD_exc_05, BR_exc_05] = bogo2_alfa7_RES_fun(BC_moded, BD, BR.a_axi, BR.a_tan); %#ok<ASGLU>
        % Set magnitude of gust and recompute forces
        BC_moded.n_sigma_exceedance = 1.0;
        [RES_axi_exc_10, RES_tan_exc_10, BD_exc_10, BR_exc_10] = bogo2_alfa7_RES_fun(BC_moded, BD, BR.a_axi, BR.a_tan); %#ok<ASGLU>
        % Set magnitude of gust and recompute forces
        BC_moded.n_sigma_exceedance = 1.5;
        [RES_axi_exc_15, RES_tan_exc_15, BD_exc_15, BR_exc_15] = bogo2_alfa7_RES_fun(BC_moded, BD, BR.a_axi, BR.a_tan); %#ok<ASGLU>
        % Set magnitude of gust and recompute forces
        BC_moded.n_sigma_exceedance = 1.9;
        [RES_axi_exc_20, RES_tan_exc_20, BD_exc_20, BR_exc_20] = bogo2_alfa7_RES_fun(BC_moded, BD, BR.a_axi, BR.a_tan); %#ok<ASGLU>
        % Extract CQ_RBM in gusts
        CQ_axi_RBM_00  = BR_exc_00.CQ_axi_root_bending_moment;
        CQ_axi_RBM_05  = BR_exc_05.CQ_axi_root_bending_moment;
        CQ_axi_RBM_10  = BR_exc_10.CQ_axi_root_bending_moment;
        CQ_axi_RBM_15  = BR_exc_15.CQ_axi_root_bending_moment;
        CQ_axi_RBM_20  = BR_exc_20.CQ_axi_root_bending_moment;
        CQ_axi_RBM_exp = BR.CQ_axi_root_bending_moment;
        % Aggregate in a vector
        CQ_axi_RBM_vector  = [CQ_axi_RBM_00, CQ_axi_RBM_05, CQ_axi_RBM_10, CQ_axi_RBM_15, CQ_axi_RBM_20, CQ_axi_RBM_exp];
        % Filter out any non-converged values
        CQ_axi_RBM_vector  = CQ_axi_RBM_vector(isfinite(CQ_axi_RBM_vector));
        % Extract CQ_RBM in gust
        CQ_axi_RBM  = max(CQ_axi_RBM_vector);
    end
end

% % 
function [BC_moded] = modify_BEMcase(BC_base, A_chord, A_theta)    
    % Make basis of 8th degree Bernstein Polynomials
    cst18 = @(mu) nchoosek(7,0) * mu.^(0) .* (1-mu).^(7);
    cst28 = @(mu) nchoosek(7,1) * mu.^(1) .* (1-mu).^(6);
    cst38 = @(mu) nchoosek(7,2) * mu.^(2) .* (1-mu).^(5);
    cst48 = @(mu) nchoosek(7,3) * mu.^(3) .* (1-mu).^(4);
    cst58 = @(mu) nchoosek(7,4) * mu.^(4) .* (1-mu).^(3);
    cst68 = @(mu) nchoosek(7,5) * mu.^(5) .* (1-mu).^(2);
    cst78 = @(mu) nchoosek(7,6) * mu.^(6) .* (1-mu).^(1);
    cst88 = @(mu) nchoosek(7,7) * mu.^(7) .* (1-mu).^(0);

    % Make shape function for twist distribution
    S_theta = @(mu) A_theta(1)*cst18(mu) + ...
                    A_theta(2)*cst28(mu) + ...
                    A_theta(3)*cst38(mu) + ...
                    A_theta(4)*cst48(mu) + ...
                    A_theta(5)*cst58(mu) + ...
                    A_theta(6)*cst68(mu) + ...
                    A_theta(7)*cst78(mu) + ...
                    A_theta(8)*cst88(mu) ;
                
    S_chord = @(mu) A_chord(1)*cst18(mu) + ...
                    A_chord(2)*cst28(mu) + ...
                    A_chord(3)*cst38(mu) + ...
                    A_chord(4)*cst48(mu) + ...
                    A_chord(5)*cst58(mu) + ...
                    A_chord(6)*cst68(mu) + ...
                    A_chord(7)*cst78(mu) + ...
                    A_chord(8)*cst88(mu) ;
    
    % Copy base BEM case to create modified BEM case
    BC_moded = BC_base;
    % Modify it
    theta_class_function_offset = 12; % gives minimum play of +-6 deg with [lb, ub] = [0.5, 1.5]
    BC_moded.theta_deg_of_mu_fun = @(mu) (BC_base.theta_deg_of_mu_fun(mu) + theta_class_function_offset).* S_theta(mu) - theta_class_function_offset;
    BC_moded.c_over_R_of_mu_fun  = @(mu)  BC_base.c_over_R_of_mu_fun( mu)                               .* S_chord(mu);
end

