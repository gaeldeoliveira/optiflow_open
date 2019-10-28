% Script to study the accumulation of a passive scalar from the convection
% of a forcing field
%
% This is the equation for the transport of a passive scalar:
%   dcdt + (u*ddx + v*ddy)c = D*lap(c)
% Decompose c into c = u_bar + u_tilde, so that:
%   ddt(u_bar+u_tilde) + (u*ddx + v*ddy)(u_bar+u_tilde) = D*lap(u_bar+u_tilde)
% Where u_tilde = 0 for all x,y at t=0 and u_bar(x,y,t) is prescribed in
% space and time. Reworking the equation, leads to:
%   ddt(u_tilde)= - ddt(u_bar) - (u*ddx + v*ddy)(u_bar+u_tilde)
% Which is a field (x,y )of ODEs, and can therefore be integrated in time.
% This is what we will do  with a simple scheme in this  script.
% The equation can be expanded with the distributive property
%   ddt(u_tilde)= - ddt(u_bar) - (u*ddx + v*ddy)(u_bar+u_tilde) + D*lap(u_bar+u_tilde)
%   ddt(u_tilde)= - ddt(u_bar) - (u*ddx + v*ddy)(u_bar) - (u*ddx + v*ddy)(u_tilde) + D*lap(u_bar) + D*lap(u_tilde)
%
% Here are our variables:
%   The convecting field is provided by the induction of a vortex system:
%   u       is the x component of the convecting field
%   v       is the y component of the convecting field
%   The forcing field u_bar is constant in time, and has a simple form for
%   development purposes: 
%   u_bar   = @(x, y, t) alpha * y
%   Where alpha is a prescribed constant (slope). In the future u_bar will
%   be provided by the Swafford profile, and its parameters will be allowed
%   to vary in time.
%   Furthermore, in later versions, dt will be related to u_bar and allowed
%   to vary in space, but for now, let us simply fix it.
%
%   Same as alfa3, except that now the time step is not constant anyu
%
%   set(findall(gcf,'-property','FontSize'),'FontSize',18)
%



%% Preparation
clear all; close all;clc; 

% Add Paths
addpath closure_relations/
addpath definitions/

%% Inputs

% Set Outer Flow Parameters
rho_inf     = 1.225;                  % [kg/m3]   Free Stream Density
nu_inf      = 1.461e-5;               % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
u_inf       = 15.16;                  % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)
msq_inf     = 0;                      % [adim.]   Unperturbed Free Stream Mach Number
u_inf_over_nu_inf = u_inf / nu_inf;   % [1/m]     Reynolds Ratio (Re number per unit lenght!)

% Set Shear Flow Parameters
hk0       = 1.35;                      % [adim.] Shape Factor
rt0       = 2466;                     % [adim.] Reynolds Theta of Starting Flow
delta0    = ... % 25e-3; 
 delta_99(hk0,rt0,u_inf_over_nu_inf); % [m    ] Boundary Layer Thickness

% Set Vortex Generator Parameters
D_VG      =  30.0e-3;                 % [m    ] distance between center of vortex pairs
d_VG      =  12.5e-3;                 % [m    ] distance between VGs in a single pair
c_VG      =  12.5e-3;                 % [m    ] chord of VG  (5cm)
h_VG      =   5.0e-3;                 % [m    ] height of VG
AoA_VG    = -18     ;                 % [deg  ] VG angle of attack (geometric, to free-stream)

% Define integration parameters
D_u_bar   = 0.000;                    % Artificial Viscosity on Forcing Field (set to 0 to avoid very high diffusion near wall, unstabilizing)
D_u_tilde = 0.001;                    % Artificial Viscosity on C_tilde Field (0.002 seems nice for 100*100)
x0        = 0;                        % Starting Time
x_end     = 0.3;%;0.5;                      % End Time
dx        = 0.0000125;                % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
dx_info   = 0.0250000;                % Time between information events
information_period=round(dx_info/dx); % Steps between information events
user_prompt = false;

% Set Outer Flow Distribution Parameters
xi_vector  = linspace(x0,x_end,100);             % Reinterpolated, number of points (100) not so relevant
ue_vector  = u_inf   * ones(size(xi_vector));   % Set flat pressure gradient with u_inf outer speed
msq_vector = msq_inf * ones(size(xi_vector));   % Set all flow as purely incompressible


%% Initialization: Create Boundary Layer Object
BL = boundary_layer6();

% Set nue 
BL.u_inf_over_nu_inf = u_inf_over_nu_inf;

% Set forcing terms (order matters! must be done after u_inf_over_nu_inf was set)
BL.set_forcing_terms(xi_vector, ue_vector, msq_vector);

%% Initialization: Create Vortex Generator Object

VG = vortex_generator(d_VG, c_VG, h_VG, AoA_VG);


%% Initialization: Create Vortex System
n_cells = 5;
S       = D_VG/2;

VD = vortex_descriptor(S, n_cells);
VD.induction_function = 3;              % Try Lamb Vortex

%% Initialization: Define forcing field

% Define forcing field
% (linear profile up to y=delta, constant 1 above, independent of x and t)
% (0*x is only there for consistent vectorization)
u_bar_fun     = @(x, y, t) ((y/delta0).^exp_BL) .*(y <  delta0) + ...
                              ones(size(y))     .*(y >= delta0) + ...
                              0*x; 

%% Make z-y Meshes

CM = crossflow_mesh(S, delta0);

z_mesh = CM.z_mesh;
y_mesh = CM.y_mesh;


%% Make Pure Shear Field
SF = shear_field(CM);


%% Make Mixed Field

% Create mixed_field (u_tilde) object 
MF = mixed_field(CM);


%% Define Initial Conditions for Integration

% Set initial position/time as current
state.x  = x0;

% % Generate Initial Boundary Layer State
% Make initial dstr0 and theta0 from initial hk0 and rt0
theta0 = rt0 / (u_inf/nu_inf);
% Make initial guess for shear stress coefficient ctau
BL.set_equilibrium_initial_conditions(theta0, hk0);          % Initial Shear Stress Coefficient determined as if BL was in equilibrium (which more less never true!!!)
ctau0 = BL.ctau0;
% Primary Integral BL Variables
state.theta = theta0;
state.hk = hk0;
state.ctau = ctau0 * 1.5;
% Secondary (dependent) Integral BL variables
state.ue        = BL.ue_function(x0);                       % Edge velocity (the BL object handles reinterpolation of the forcing field, for historical reasons!)
state.rt        = state.theta * (state.ue / nu_inf);        % Should be equal to rt0


% % Generate Initial State of Vortex Descriptors
% Differs with the type of chosen vortex core
if     VD.induction_function == 2
    % Rankine Vortex Induction
    [y_v0, z_v0, gamma_v0, sigma_v0] = VG.initial_strenght_stripline_model(u_inf);
elseif VD.induction_function == 3
    % Lamb Vortex Induction
    % Use classical prandtl for gamma
    %[   ~,    ~,  gamma_v0, ~] = VG.initial_strenght_prandtl_model(u_inf);
    % Or smart alternative
    u_bar_v0    = SF.u_bar_at_y(state.hk, state.rt, state.ue, nu_inf, VG.h_VG);
    [   ~,    ~,  gamma_v0, ~] = VG.initial_strenght_prandtl_model(u_bar_v0);
    % And Wendt2001 for peak vorticity (something is wrong with our implementation of the gamma)
    [y_v0, z_v0, ~ , sigma_v0] = VG.initial_strenght_wendt2001_model(u_inf, delta_99(hk0,rt0,u_inf_over_nu_inf));
end
% Store into state variables
state.y_v     = y_v0;             % [m /s ] distance of filament core center to wall
state.z_v     = z_v0;             % [m /s ] distance of filament core center to center of Vg pair
state.gamma_v = gamma_v0;         % [m /s ] circulation of vortex filament
state.sigma_v = sigma_v0;         % [m    ] core size of vortex filament (Rankine)
                                  % [1 /s ] peak vorticity of vortex filament (Lamb)
% % Generate Initial Age of Vortex
nu_v          = nu_inf;
state.t_v     = initial_vortex_wendt1995_model(VG, gamma_v0, sigma_v0, nu_v);

% Generate initial u_tilde mesh (set at 0 for initial state)
u_tilde0 = MF.initial_u_tilde();
state.u_tilde  = u_tilde0;

% Make a copy of the original state for future reference
state0 = state;



%% Create Field Mixer Object
FM = field_mixer(VD, CM, SF, MF, D_u_bar, D_u_tilde);

%% Test thicknesses integration
SMI = shear_mixed_integrator(CM, SF, MF);
[dstr_bar, dstr_tilde, dstr_hat] = SMI.dstr_direct(state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
[theta_bar, theta_tilde, theta_hat] = SMI.theta_direct(state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
[delta3_bar, delta3_tilde, delta3_hat] = SMI.delta3_direct(state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
h_bar = dstr_bar / theta_bar;
h_hat = dstr_hat / theta_hat;
disp(['h_bar = ' num2str(h_bar) '    h_hat = ' num2str(h_hat)]);

%% Prepare Plotting
% Open Figures
figure(1) % First   plot (Mixed Field)
figure(2) % Second  plot (Mixing Field)
figure(3) % Third   plot (Composite Plot)
figure(4) % Fourth  plot (Velocity Profiles)
figure(5) % Fifth   plot (Core Position)
figure(6) % Sixth   plot (Integral BL Parameters)
figure(7) % Seventh plot (Core Diffusion)
figure(8) % Eight   plot (Composite Plot)

%% Prepare File Output
% Make Output Filena and Folder Strings 
painter_code  = '-dpdf';
painter_code2 = '-dpng';
fig_folder    = [ datestr(now, 30), '/figures/'];
mat_folder    = [ datestr(now, 30), '/matfiles/'];
file_ext      = '.pdf';
file_ext2     = '.png';
painter_opt2  = '-r300';

% Make folders
mkdir(mat_folder);
mkdir(fig_folder);


%% Now Load Experimental data for Comparison
% Define Filename Root
fname_root     = 'ZPG_yaw0_VGrect_d';

% Define Streamwise Stances for Reading
d_list         = [ 5       ,  6       ,  7       ,  8       ,  9       , 10       , 25       , 50       ];
% Mark Stances for Filtering
d_list_flag    = [ 1       ,  1       ,  1       ,  1       ,  1       ,  1       ,  0       ,  0       ];
% Define Filters
z_min_list     = [-0.007549, -0.007454, -0.007171, -0.007454, -0.007738, -0.007738];
z_max_list     = [ 0.00944 ,  0.009157,  0.009723,  0.009629,  0.009912,  0.009912];
y_min_list     = [ 0.0134  ,  0.01708 ,  0.02152 ,  0.0252  ,  0.02973 ,  0.02973 ];
y_max_list     = [ 0.02756 ,  0.03152 ,  0.03388 ,  0.03388 ,  0.03388 ,  0.03388 ];
% Read Experimental Data
ED             = BLT_experiment_data_reader(fname_root, d_list, h_VG);
% Set filters
ED.set_filters(d_list_flag, z_min_list, z_max_list, y_min_list, y_max_list);
% Filter data
ED.filter_data();

% % Set Scalling 
% Set Crossflow scaling variables
ED.set_crossflow_scaling_variables(S, delta0);


%% Run Evolution Loop

% Now start cycle
for n = 0:round(x_end/dx)
    
    % % Evolve Boundary Layer (Drela 1987)
    % Get spatial rate of change of primary variables
    [d_theta_dx, d_h_dx, d_ctau_dx] = BL.manual_fun_wrapper(state.theta, state.hk, state.ctau, state.x);
    % March parameters (primary BL)
    new_state.theta     = state.theta + d_theta_dx * dx;
    new_state.hk        = state.hk    + d_h_dx     * dx;
    new_state.ctau      = state.ctau  + d_ctau_dx  * dx;
    
    % % Filament Position Evolution (Jones 1957)
    % Determine induction at (primary) vortex core center
    [w_v, v_v] = VD.induction_vortex_center(state.z_v, state.y_v, state.gamma_v, state.sigma_v);
    % Determine convection speed at (primary) vortex core center
    u_bar_v     = SF.u_bar_at_y(state.hk, state.rt, state.ue, nu_inf, state.y_v);
    % Determine lagragiant time derivative to space (x) at core center
    dT_v_dx = 1 / u_bar_v;        % don't confuse with dT_dx of u_tilde field! 
                                 % same stuff, but for different points (here single point, there mesh)
    % Determine spatial rate of change of filament position
    dy_v_dx       = v_v * dT_v_dx;
    dz_v_dx       = w_v * dT_v_dx; 
    % March filament position in space
    new_state.y_v = state.y_v + dy_v_dx * dx;
    new_state.z_v = state.z_v + dz_v_dx * dx;
    
    % % Vortex Core Diffusion Evolution (Squirre 1965)
    % Determine rate of diffusion in time
    dsigma_v_dt_v = - (1 / (4*pi())) * (state.gamma_v ./ nu_v) * (1 ./ state.t_v).^2;
    % Determine convection speed at (primary) vortex core center
    u_bar_v     = SF.u_bar_at_y(state.hk, state.rt, state.ue, nu_inf, state.y_v);
    % Determine lagragiant time derivative to space (x) at core center
    dT_v_dx = 1 / u_bar_v;        % don't confuse with dT_dx of u_tilde field! 
                                  % same stuff, but for different points (here single point, there mesh)
    % Determine rate of diffusion in space
    dsigma_v_dx = dsigma_v_dt_v * dT_v_dx;
    % March Core Diffusion Equations
    new_state.t_v     = state.t_v     + dT_v_dx     * dx;
    new_state.sigma_v = state.sigma_v + dsigma_v_dx * dx;
    
%   Squirre single unit test!
%     dsigma_v_dt_v = - (1 / (4*pi())) * (gamma_v ./ nu_v) * (1 ./ t_v).^2;
%     u_bar_v     = SF.u_bar_at_y(state.hk, state.rt, state.ue, nu_inf, state.y_v);
%     dT_v_dx = 1 / u_bar_v;        
%     dsigma_v_dx = dsigma_v_dt_v * dT_v_dx;
%     t_v     = t_v     + dT_v_dx     * dx
%     sigma_v = sigma_v + dsigma_v_dx * dx
    
    
    % % Mixing Field Evolution (Oliveira 2016)
    % Filter and Compute spatial rate of change of u_tilde
    [state.u_tilde, u_bar, du_tilde_dx, dT_dx, mag_mesh] = ...
                    FM.u_tilde_rate_of_change( state, nu_inf);
    % March u_tilde in space
    new_state.u_tilde = state.u_tilde + du_tilde_dx  * dx;
    
    % % March stuff, Transmit static States, Make it Consistent
    % Position
    new_state.x       = state.x + dx;
    
    % Update dependent parameters (secondary BL)
    new_state.ue        = BL.ue_function(new_state.x);
    new_state.rt        = new_state.theta * (new_state.ue / nu_inf);
    
    % Transmit (currently) static parameters to new state
    new_state.gamma_v   = state.gamma_v;
    % new_state.sigma_v   = state.sigma_v;
    
    % Update state into new_state! 
    % state = new_state; Move this after plotting
    
    % Information Output
    if (n/ information_period == round(n / information_period))        
        % CFL Estimation
        c = max(max(mag_mesh));
        dh_min = min(min(diff(CM.z_range)) , min(diff(CM.y_range)));
        dh_max = max(max(diff(CM.z_range)) , max(diff(CM.y_range)));
        dt_max = max(max(dT_dx*dx));
        CFL = c*dt_max / dh_min;    
        PE  = c*dh_max / D_u_tilde; % U *deltaX / D
        disp(['n = ' num2str(n) '    x = ' num2str(state.x) '    CFL = ' num2str(CFL), '    PE = ' num2str(PE)])
        
        
        % % % Plotting and File Output
        % % Update Output Identifier String        
        id_str = ['_' num2str(round(state.x / h_VG*10)/10) 'h'];
        h_str  = [num2str(round(state.x / h_VG*10)/10) ' h^{VG}'];
        
        % % Plot Mixed Flow Velocity Field (u_tilde)
        set(0, 'CurrentFigure', 1)
        surf(z_mesh/S, y_mesh/delta0, state.u_tilde / state.ue); hold on;
        offset = 1; contour3(z_mesh/S, y_mesh/delta0, state.u_tilde / state.ue, [-0.1 -0.05 0.05 0.1] + offset ,'w');
        xlabel('z/S'); ylabel('y/\delta'); colorbar;
        title(['u^{tilde} / u_e    at x = ' h_str]);
        view(2); shading flat;  caxis([-0.15 0.15]);
        axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        % print(painter_code , [fig_folder, 'mixed_flow', id_str , file_ext ]);
        print(painter_code2, [fig_folder, 'mixed_flow', id_str , file_ext2], painter_opt2);
        
        % % Plot Vortical Flow Component (|v_tilde+w_tilde|)
        set(0, 'CurrentFigure', 2)
        [~, ~, mag_mesh] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, 1);
        surf(z_mesh/S, y_mesh/delta0, mag_mesh / new_state.ue); hold on;
        offset = 1; contour3(z_mesh/S, y_mesh/delta0, mag_mesh / new_state.ue + offset, [0.05 0.1 0.2 0.4] + offset ,'w');
        xlabel('z/S'); ylabel('y/\delta'); colorbar;
        title(['|v^{tilde} + w^{tilde}|/ u_e   at x = ' h_str]);
        view(2); shading flat;
        axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
        caxis([0 1]); hold off;
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        %print(painter_code , [fig_folder, 'vortical_flow', id_str , file_ext ]);
        print(painter_code2, [fig_folder, 'vortical_flow', id_str , file_ext2], painter_opt2);
             
        % % Composite Plot of Mixed and Vortical Fields
        set(0, 'CurrentFigure', 3)
        % Mixed Flow Velocity Field (u_tilde)
        subplot(221)
        surf(z_mesh/S, y_mesh/delta0, state.u_tilde / state.ue); hold on;
        offset = 1; contour3(z_mesh/S, y_mesh/delta0, state.u_tilde / state.ue, [-0.1 -0.05 0.05 0.1] + offset ,'w');
        xlabel('z/S'); ylabel('y/\delta'); colorbar;
        title(['u^{tilde} / u_e    at x = ' h_str]);
        view(2); shading flat; caxis([-0.15 0.15]);
        axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);               
        % Vortical Flow Field (|v_tilde+w_tilde|)
        subplot(222)
        [w_mesh, v_mesh, mag_mesh] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, 1);
        surf(z_mesh/S, y_mesh/delta0, mag_mesh / new_state.ue); hold on
        offset = 1; contour3(z_mesh/S, y_mesh/delta0, mag_mesh / new_state.ue + offset, [0.05 0.1 0.2 0.4] + offset ,'w');
        xlabel('z/S'); ylabel('y/\delta'); colorbar;
        title(['|v^{tilde} + w^{tilde}|/ u_e   at x = ' h_str]);
        view(2); shading flat; 
        axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
        caxis([0 1]); hold off;        
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        % print(painter_code , [fig_folder, 'composite_mixed_vortical', id_str , file_ext ]);
        print(painter_code2, [fig_folder, 'composite_mixed_vortical', id_str , file_ext2], painter_opt2);
                               
        % % Composite Plot of Mixed and Total Fields
        set(0, 'CurrentFigure', 8)
        % Total Flow Field (u_tilde+u_bar|)
        subplot(222)
        u_bar_mesh = SF.u_bar_over_mesh(state.hk, state.rt, state.ue, state.ue / nu_inf);
        surf(z_mesh/S, y_mesh/delta0, (state.u_tilde + u_bar_mesh )/ state.ue)
        xlabel('z/S'); ylabel('y/\delta'); colorbar;
        title(['Predicted (u^{tilde}+u^{bar})/u_e   at x = ' h_str]);
        view(2); shading flat; caxis([0 1]);
        axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);  
        % Total Flow Field (Experimental)
        subplot(221)
        [z_exp_mesh, y_exp_mesh, u_exp_mesh] = ED.interp_u_mesh(state.x/h_VG);
        surf(z_exp_mesh/S, y_exp_mesh/delta0, u_exp_mesh/ state.ue)
        xlabel('z/S'); ylabel('y/\delta'); colorbar;
        title(['Experimental u/u_e   at x = ' h_str]);
        view(2); shading flat; caxis([0 1]);
        axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);  
        
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        % print(painter_code , [fig_folder, 'composite_mixed_total', id_str , file_ext ]);
        print(painter_code2, [fig_folder, 'composite_mixed_total', id_str , file_ext2], painter_opt2);
        
        % % Plot Velocity Profiles
        % Centerline
        set(0, 'CurrentFigure', 4)
        subplot(121)
        j_center = CM.j_center;
        plot( new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
                          u_bar(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
                          u_bar(:, j_center )/ new_state.ue + new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
                          ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, 0) / u_inf, CM.y_range/delta0, 'k');
        grid on; axis([-0.5 1.5 0 2]);
        title('Central Symmetry Line');
        xlabel('u/u_{inf}'); ylabel('y/\delta');
        % Mid Line
        subplot(122)
        j_quarter = round(j_center/2);
        plot(new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
                         u_bar(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
                         u_bar(:, j_quarter)/ new_state.ue + new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
                         ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, -1) / u_inf, CM.y_range/delta0, 'k');
        grid on; axis([-0.5 1.5 0 2]);
        title('Side Symmetry Line');
        xlabel('u/u_{inf}'); ylabel('y/\delta');
        legend('Mixed Flow' , 'Pure Shear Flow', 'Total Predicted Flow', 'Experimental Result')
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        print(painter_code , [fig_folder, 'velocity_profiles', id_str , file_ext ]);
        print(painter_code2, [fig_folder, 'velocity_profiles', id_str , file_ext2], painter_opt2);
        
        
        % % Replot Figure 4 on smaller scale
        set(0, 'CurrentFigure', 5)
        subplot(221)
        j_center = CM.j_center;
        plot( new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
                          u_bar(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
                          u_bar(:, j_center )/ new_state.ue + new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
                          ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, 0) / u_inf, CM.y_range/delta0, 'k');
        grid on; axis([-0.5 1.5 0 1]);
        title(['Central Sym. Line at x = ' h_str]);
        xlabel('u/u_{inf}'); ylabel('y/\delta');
        % Mid Line
        subplot(222)
        j_quarter = round(j_center/2);
        plot(new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
                         u_bar(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
                         u_bar(:, j_quarter)/ new_state.ue + new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
                         ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, -1) / u_inf, CM.y_range/delta0, 'k');
        grid on; axis([-0.5 1.5 0 1]);
        title(['Side Sym. Line at x = ' h_str]);
        xlabel('u/u_{inf}'); ylabel('y/\delta');
        legend('Mixed Flow' , 'Pure Shear Flow', 'Total Predicted Flow', 'Experimental Result')
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        print(painter_code , [fig_folder, 'velocity_profiles_subplot', id_str , file_ext ]);
        print(painter_code2, [fig_folder, 'velocity_profiles_subplot', id_str , file_ext2], painter_opt2);
        
        
        
%         % % Plot Vortex Core Center Positions (Generalized Jones (1957) model)
%         set(0, 'CurrentFigure', 5)
%         %plot3( state.x ,  state.z_v/S    , state.y_v/delta0, 'ok'); hold on; grid on;
%         %plot3( state.x , -state.z_v/S    , state.y_v/delta0, 'or');
%         %plot3( state.x ,  state.z_v/S - 2, state.y_v/delta0, 'ok');
%         %plot3( state.x , -state.z_v/S + 2, state.y_v/delta0, 'or');       
%         %xlabel('x/L')  ;  ylabel('z/S')  ; zlabel('y/\delta');
%         %axis([x0-eps state.x min(CM.z_range)/S max(CM.z_range)/S min(CM.y_range)/delta0 max(CM.y_range)/delta0])
%         %title('Vortex Core Center')
%         subplot(221)
%         % Plot balls, color indicates pair of belonging
%         plot(  state.z_v/S    , state.y_v/delta0, 'ok'); hold on; grid on;
%         plot( -state.z_v/S    , state.y_v/delta0, 'ok');
%         plot(  state.z_v/S - 2, state.y_v/delta0, 'or');
%         plot( -state.z_v/S + 2, state.y_v/delta0, 'or');       
%         % Plot symbols indicating orientation, color matches belonging ball
%         cw_symbol   = plot(  state.z_v/S    , state.y_v/delta0, 'xk'); hold on; grid on;
%         ccw_symbol  = plot( -state.z_v/S    , state.y_v/delta0, '+k');
%                      plot(  state.z_v/S - 2, state.y_v/delta0, 'xr');
%                      plot( -state.z_v/S + 2, state.y_v/delta0, '+r');       
%         sym_symbol  = plot( [0 0], [min(CM.y_range)/delta0 max(CM.y_range)/delta0], 'k-.');
%         side_symbol = plot( [1 1], [min(CM.y_range)/delta0 max(CM.y_range)/delta0], 'r-.');
%                       plot(-[1 1], [min(CM.y_range)/delta0 max(CM.y_range)/delta0], 'r-.');
%         xlabel('z/S')  ; ylabel('y/\delta');
%         axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
%         legend([cw_symbol, ccw_symbol, sym_symbol, side_symbol], '\Gamma<0', '\Gamma>0', 'Pair', 'Side');
%         legend([cw_symbol, ccw_symbol], '\Gamma<0', '\Gamma>0');
%         subplot(222)
%         plot( state.x / h_VG , state.y_v/delta0, 'ok'); hold on; grid on;
%         xlabel('x/h^{VG}')  ;  ylabel('y/\delta');
%         % Output to File
%         set(gcf , 'PaperType', 'A5');
%         orient landscape;
%         print(painter_code , [fig_folder, 'core_center', id_str , file_ext ]);
%         print(painter_code2, [fig_folder, 'core_center', id_str , file_ext2], painter_opt2);
        
              
        % % Document closures from shear_mixed_integrator
        [  dstr_bar,   dstr_tilde,   dstr_hat] = SMI.dstr_direct(  state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
        [ theta_bar,  theta_tilde,  theta_hat] = SMI.theta_direct( state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
        [delta3_bar, delta3_tilde, delta3_hat] = SMI.delta3_direct(state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
        h_bar  =   dstr_bar / theta_bar;  h_hat =   dstr_hat / theta_hat;
        hs_bar = delta3_bar / theta_bar; hs_hat = delta3_hat / theta_hat;
        disp(['     h_bar = ' num2str(    h_bar)  '         h_hat = ' num2str(     h_hat)]);
        disp(['    hs_bar = ' num2str(   hs_bar)  '        hs_hat = ' num2str(    hs_hat)]);
        disp(['  dstr_bar = ' num2str( dstr_bar)  '      dstr_hat = ' num2str(  dstr_hat) '      dstr_tilde = ' num2str(  dstr_tilde)]);
        disp([' theta_bar = ' num2str( theta_bar) '     theta_hat = ' num2str( theta_hat) '     theta_tilde = ' num2str( theta_tilde)]);
        disp(['delta3_bar = ' num2str(delta3_bar) '    delta3_hat = ' num2str(delta3_hat) '    delta3_tilde = ' num2str(delta3_tilde)]);
        
        % % Plot Integral Boundary Layer Parameters (Pure Shear and Total Flow)
        set(0, 'CurrentFigure', 6)
        % Shape Factor
        subplot(211)
        plot( state.x/h_VG, h_bar   , 'ok'); hold on; grid on;
        plot( state.x/h_VG, h_hat   , 'or');
        xlabel('x/h^{VG}'); ylabel('H_k - Shape factor')
        title('Shape Factor')
        legend('Pure Shear Flow', 'Total Flow')
        axis([x0-eps state.x/h_VG 1.3 1.5]);
        % Displacement Thickness
        subplot(212)
        plot( state.x / h_VG, dstr_bar   , 'ok'); hold on; grid on;
        plot( state.x / h_VG, dstr_hat   , 'or');
        plot( state.x / h_VG, dstr_tilde , 'ob');
        legend('Pure Shear Flow' , 'Total Flow', 'Mixed Flow')
        xlabel('x/h^{VG}'); ylabel('\delta* - Displacement Thickness')
        title('Displacement Thickness')
        legend('Pure Shear Flow' , ' Total Flow', 'Mixed Flow')
        axis([x0-eps state.x 0 2*(hk0*rt0/(u_inf/nu_inf))]);
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        print(painter_code , [fig_folder, 'IBL_parameters', id_str , file_ext ]);
        print(painter_code2, [fig_folder, 'IBL_parameters', id_str , file_ext2], painter_opt2);
        
        % % Plot Circulation and Peak Vorticity (Squirre (1965) Model)
        % Peak Vorticity
        set(0, 'CurrentFigure', 7)
        subplot(211)
        plot( state.x / h_VG, state.sigma_v , 'ok'); hold on; grid on;
        xlabel('x/h^{VG}'); ylabel('\omega^{max}_{0}')
        title('Peak Vorticity in Lamb Vortex Core')
        axis([x0-eps state.x/h_VG 0 state0.sigma_v*1.2]);
        % Circulation
        subplot(212)
        plot( state.x / h_VG, state.gamma_v , 'ok'); hold on; grid on;
        xlabel('x/h^{VG}'); ylabel('\Gamma_{v}')
        title('Filament Circulation per Unit Lenght')
        axis([x0-eps state.x/h_VG 0 state0.gamma_v*1.2]);
        % Output to File
        set(gcf , 'PaperType', 'A5');
        orient landscape;
        print(painter_code , [fig_folder, 'core_diffusion', id_str , file_ext ]);
        print( gcf, painter_code2, [fig_folder, 'core_diffusion', id_str , file_ext2], painter_opt2);
        
        % Save Current Workspace State
        save([mat_folder, 'state' , id_str, '.mat']);
        
        
        % Update Graphical Output
        drawnow();              
        % Prompt for continuation
        if(user_prompt==true)
            input('Go on?');
        else
            close(1);
            close(2);
            close(3);
            close(8);
            drawnow();
            figure(1);
            figure(2);
            figure(3);
            figure(8);
        end
    end
    
    % Update state into new_state!
    state = new_state;
end


% Later on, remember to export some images
% Print it
% set(gcf , 'PaperType', 'A5');
% orient landscape;
% print ('-dpdf', ['wendt_cases_subplot.pdf']) %#ok<NBRAK>




