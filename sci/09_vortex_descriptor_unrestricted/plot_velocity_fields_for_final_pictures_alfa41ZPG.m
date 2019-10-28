%% Preparation
clear all; close all;clc; 

% Add Paths
addpath closure_relations/
addpath definitions/
addpath intersections/

%% Inputs

% Set Outer Flow Parameters
rho_inf     = 1.225;                  % [kg/m3]   Free Stream Density
nu_inf      = 1.461e-5;               % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
u_inf       = 15.16;                  % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)
msq_inf     = 0;                      % [adim.]   Unperturbed Free Stream Mach Number
u_inf_over_nu_inf = u_inf / nu_inf;   % [1/m]     Reynolds Ratio (Re number per unit lenght!)

% Set Shear Flow Parameters
hk0       = 1.41;                     % [adim.] Shape Factor
rt0       = 2498;                     % [adim.] Reynolds Theta of Starting Flow
delta0    = ... % 25e-3; 
 delta_99(hk0,rt0,u_inf_over_nu_inf); % [m    ] Boundary Layer Thickness

% Set Vortex Generator Parameters
D_VG      =  30.0e-3;                 % [m    ] distance between center of vortex pairs
d_VG      =  12.5e-3;                 % [m    ] distance between VGs in a single pair
c_VG      =  12.5e-3;                 % [m    ] chord of VG  (5cm)
h_VG      =   5.0e-3;                 % [m    ] height of VG
AoA_VG    = -18     ;                 % [deg  ] VG angle of attack (geometric, to free-stream)

% Define integration parameters
D_u_bar   =    0.000;                 % Artificial Viscosity on Forcing Field (set to 0 to avoid very high diffusion near wall, unstabilizing)
D_u_tilde =    0.002;                 % Artificial Viscosity on C_tilde Field (0.002 seems nice for 100*100)
x0        =    -0.50 * h_VG;          % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
x_end     =    51    * h_VG;%0.25;%;0.5;        % End Time (0.25 = 50 * h_vg)
dx        =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
dx_info   =    0.50 * h_VG;           % Time between information events
information_period=round(dx_info/dx); % Steps between information events

% Define output parameters
user_prompt = false;
img_output = true;
mat_output = true;
baldwin_factor = 2;

% Set new step size
pod_dx                 = 1.25e-3;     % Up to 1000 times bigger than the dx of the normal stuff!
% Set pod information period
pod_information_period = round(dx_info/pod_dx);


% Set Outer Flow Distribution Parameters
xi_vector  = linspace(x0,x_end,100);            % Reinterpolated, number of points (100) not so relevant
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
VD.induction_function = 3;              % Try Lamb Vortex (had something wrong with damping, was falling back to singular vortex!)
VD.induction_function = 4;              % Try Debugged Lamb Vortex

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
elseif VD.induction_function == 4
    % Lamb Vortex Induction
    % Use classical prandtl for gamma
    %[   ~,    ~,  gamma_v0, ~] = VG.initial_strenght_prandtl_model(u_inf);
    % Or smart alternative
    u_bar_v0    = SF.u_bar_at_y(state.hk, state.rt, state.ue, nu_inf, VG.h_VG);
    [   ~,    ~,  gamma_v0, ~] = VG.initial_strenght_prandtl_model(u_bar_v0);
    % And Wendt2001 for peak vorticity (something is wrong with our implementation of the gamma)
    [y_v0, z_v0, gamme_v0W , sigma_v0] = VG.initial_strenght_wendt2001_model(u_inf, delta_99(hk0,rt0,u_inf_over_nu_inf));
end
% Store into state variables
state.y_v     = h_VG ;            % [m /s ] distance of filament core center to wall              % exploration = h_VG * 2/3;
state.z_v     = z_v0 ;            % [m /s ] distance of filament core center to center of Vg pair % exploration = z_v0 - 1/4 * c_VG * sin(abs(AoA_VG)* pi/180);
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
FM.baldwin_factor = baldwin_factor;

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
figure(9) % Ninth   plot (Velocity Sidecut)

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
% ED             = BLT_experiment_data_reader(fname_root, d_list, h_VG);
% % Set filters
% ED.set_filters(d_list_flag, z_min_list, z_max_list, y_min_list, y_max_list);
% % Filter data
% ED.filter_data();
% % % Set Scalling 
% % Set Crossflow scaling variables
% ED.set_crossflow_scaling_variables(S, delta0);

% Read
ED             = BLT_experiment_data_reader(fname_root, d_list, h_VG);
% Set filters
ED.set_filters(d_list_flag, z_min_list, z_max_list, y_min_list, y_max_list);
% Filter data
ED.filter_data();
% Set Crossflow scaling variables
ED.set_crossflow_scaling_variables(S, delta0);
% Extrapolate
ED.extrapolate_u_fields(u_inf, nu_inf)
% And make default
ED.make_extrapolated_fields_default()

%% And now prepare to capture snapshots
snapshots = struct();
snapshots.N_steps      = round(round(x_end/dx)/4)+1;    % Only collect one out of 4 steps to avoid running out of memory on laptop!
snapshots.z_mesh       = CM.z_mesh;
snapshots.y_mesh       = CM.y_mesh;
snapshots.u_tilde_data = zeros(size(snapshots.z_mesh,1), size(snapshots.z_mesh,2), snapshots.N_steps);
n_snapshot             = 0;

%% This is where the iterative loop would start!
close all

%% Now load some stuff

state_05h     = load('20180812T004926/matfiles/full_state_asy_5h.mat'  ); state_05h = state_05h.state;
state_10h     = load('20180812T004926/matfiles/full_state_asy_10h.mat' ); state_10h = state_10h.state;
state_25h     = load('20180812T004926/matfiles/full_state_asy_25h.mat' ); state_25h = state_25h.state;
state_50h     = load('20180812T004926/matfiles/full_state_asy_50h.mat' ); state_50h = state_50h.state;


%% And plot velocity profiles for x=5h
% Set state for x=5h
new_state = state_05h;
state     = new_state;
% Filter and Compute spatial rate of change of u_tilde (reconstruction)
[state.u_tilde, u_bar, du_tilde_dx, dT_dx, mag_mesh] = ...
                    FM.u_tilde_rate_of_change( state, nu_inf);
% Open figure
figure(405)
% Centerline
subplot(231)
j_center = CM.j_center;
h4_231   = plot( new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue + new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, 0) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Central Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Mid Line
subplot(233)
j_quarter = round(j_center/2);
h4_233    = plot(new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue + new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, -1) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Side Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(-1)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% legend('Mixed Flow' , 'Pure Shear Flow', 'Total Predicted Flow', 'Experimental Result')
legend(h4_233([4 3 2 1]), 'Experiment' , 'Reconstruction', 'Shear flow', 'Interaction flow')
% Intermediate Line
subplot(232)
j_eighth = round((j_center+j_quarter)/2);
h4_232   = plot(new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue + new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, - CM.z_range(j_eighth) ./ CM.S) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); % grid on; axis([-0.5 1.5 0 2]);
% title('Intermediate Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0.5)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Update colors
blue   = [0 , 147, 208]/255;%[0         0.4470    0.7410];
red    = [0 ,   0,   0]    ;%[0.8500    0.3250    0.0980];
yellow = [0 ,   0,   0]    ;%[0.9290    0.6940    0.1250];
black  = [0 ,   0,   0]    ;%[0         0         0     ];
% Set them
h4_231(1).Color = black ; h4_232(1).Color = h4_231(1).Color; h4_233(1).Color = h4_231(1).Color;
h4_231(2).Color = yellow; h4_232(2).Color = h4_231(2).Color; h4_233(2).Color = h4_231(2).Color;
h4_231(3).Color = red   ; h4_232(3).Color = h4_231(3).Color; h4_233(3).Color = h4_231(3).Color;
h4_231(4).Color = blue  ; h4_232(4).Color = h4_231(4).Color; h4_233(4).Color = h4_231(4).Color;
% And set linestyles
h4_231(1).LineStyle = '-.' ; h4_232(1).LineStyle = h4_231(1).LineStyle; h4_233(1).LineStyle  = h4_231(1).LineStyle;
h4_231(2).LineStyle = ':' ; h4_232(2).LineStyle  = h4_231(2).LineStyle; h4_233(2).LineStyle  = h4_231(2).LineStyle;
h4_231(3).LineStyle = '--'; h4_232(3).LineStyle  = h4_231(3).LineStyle; h4_233(3).LineStyle  = h4_231(3).LineStyle;
h4_231(4).LineStyle = '-' ; h4_232(4).LineStyle  = h4_231(4).LineStyle; h4_233(4).LineStyle  = h4_231(4).LineStyle;

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG,1), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-deps' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
savefig(        ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h', '.fig'])

%% Total Flow Field (u_tilde+u_bar|)
figure(805)
subplot(222)
u_bar_mesh = SF.u_bar_over_mesh(state.hk, state.rt, state.ue, state.ue / nu_inf);
surf(z_mesh/S, y_mesh/delta0, (state.u_tilde + u_bar_mesh )/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Predicted U/U_e   at x = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on
% Total Flow Field (Experimental)
subplot(221)
[z_exp_mesh, y_exp_mesh, u_exp_mesh] = ED.interp_u_mesh(state.x/h_VG);
surf(z_exp_mesh/S, y_exp_mesh/delta0, u_exp_mesh/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Experimental U/U_e   at X/h = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_fields_', num2str(state.x/h_VG,1), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.png'], '-r600');
savefig(        ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.fig'])
colormap(bone(16));
print('-deps' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.png'], '-r600');

%% And plot velocity profiles again for x=10h
% Set state for x=10h
new_state = state_10h;
state     = new_state;
% Filter and Compute spatial rate of change of u_tilde (reconstruction)
[state.u_tilde, u_bar, du_tilde_dx, dT_dx, mag_mesh] = ...
                    FM.u_tilde_rate_of_change( state, nu_inf);
% Open figure
figure(410)
% Centerline
subplot(231)
j_center = CM.j_center;
h4_231   = plot( new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue + new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, 0) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Central Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Mid Line
subplot(233)
j_quarter = round(j_center/2);
h4_233    = plot(new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue + new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, -1) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Side Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(-1)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% legend('Mixed Flow' , 'Pure Shear Flow', 'Total Predicted Flow', 'Experimental Result')
legend(h4_233([4 3 2 1]), 'Experiment' , 'Reconstruction', 'Shear flow', 'Interaction flow')
% Intermediate Line
subplot(232)
j_eighth = round((j_center+j_quarter)/2);
h4_232   = plot(new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue + new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, - CM.z_range(j_eighth) ./ CM.S) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); % grid on; axis([-0.5 1.5 0 2]);
% title('Intermediate Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0.5)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Update colors
blue   = [0 , 147, 208]/255;%[0         0.4470    0.7410];
red    = [0 ,   0,   0]    ;%[0.8500    0.3250    0.0980];
yellow = [0 ,   0,   0]    ;%[0.9290    0.6940    0.1250];
black  = [0 ,   0,   0]    ;%[0         0         0     ];
% Set them
h4_231(1).Color = black ; h4_232(1).Color = h4_231(1).Color; h4_233(1).Color = h4_231(1).Color;
h4_231(2).Color = yellow; h4_232(2).Color = h4_231(2).Color; h4_233(2).Color = h4_231(2).Color;
h4_231(3).Color = red   ; h4_232(3).Color = h4_231(3).Color; h4_233(3).Color = h4_231(3).Color;
h4_231(4).Color = blue  ; h4_232(4).Color = h4_231(4).Color; h4_233(4).Color = h4_231(4).Color;
% And set linestyles
h4_231(1).LineStyle = '-.' ; h4_232(1).LineStyle = h4_231(1).LineStyle; h4_233(1).LineStyle  = h4_231(1).LineStyle;
h4_231(2).LineStyle = ':' ; h4_232(2).LineStyle  = h4_231(2).LineStyle; h4_233(2).LineStyle  = h4_231(2).LineStyle;
h4_231(3).LineStyle = '--'; h4_232(3).LineStyle  = h4_231(3).LineStyle; h4_233(3).LineStyle  = h4_231(3).LineStyle;
h4_231(4).LineStyle = '-' ; h4_232(4).LineStyle  = h4_231(4).LineStyle; h4_233(4).LineStyle  = h4_231(4).LineStyle;

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG,1), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-deps' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
savefig(        ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h', '.fig'])


%% Total Flow Field (u_tilde+u_bar|)
figure(810)
subplot(222)
u_bar_mesh = SF.u_bar_over_mesh(state.hk, state.rt, state.ue, state.ue / nu_inf);
surf(z_mesh/S, y_mesh/delta0, (state.u_tilde + u_bar_mesh )/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Predicted U/U_e   at x = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on
% Total Flow Field (Experimental)
subplot(221)
[z_exp_mesh, y_exp_mesh, u_exp_mesh] = ED.interp_u_mesh(state.x/h_VG);
surf(z_exp_mesh/S, y_exp_mesh/delta0, u_exp_mesh/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Experimental U/U_e   at X/h = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_fields_', num2str(state.x/h_VG,1), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.png'], '-r600');
savefig(        ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.fig'])
colormap(bone(16));
print('-deps' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.png'], '-r600');


%% And plot velocity profiles again for x=25h
% Set state for x=25h
new_state = state_25h;
state     = new_state;
% Filter and Compute spatial rate of change of u_tilde (reconstruction)
[state.u_tilde, u_bar, du_tilde_dx, dT_dx, mag_mesh] = ...
                    FM.u_tilde_rate_of_change( state, nu_inf);
% Open figure
figure(425)
% Centerline
subplot(231)
j_center = CM.j_center;
h4_231   = plot( new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue + new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, 0) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Central Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Mid Line
subplot(233)
j_quarter = round(j_center/2);
h4_233    = plot(new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue + new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, -1) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Side Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(-1)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% legend('Mixed Flow' , 'Pure Shear Flow', 'Total Predicted Flow', 'Experimental Result')
legend(h4_233([4 3 2 1]), 'Experiment' , 'Reconstruction', 'Shear flow', 'Interaction flow')
% Intermediate Line
subplot(232)
j_eighth = round((j_center+j_quarter)/2);
h4_232   = plot(new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue + new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, - CM.z_range(j_eighth) ./ CM.S) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); % grid on; axis([-0.5 1.5 0 2]);
% title('Intermediate Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0.5)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Update colors
blue   = [0 , 147, 208]/255;%[0         0.4470    0.7410];
red    = [0 ,   0,   0]    ;%[0.8500    0.3250    0.0980];
yellow = [0 ,   0,   0]    ;%[0.9290    0.6940    0.1250];
black  = [0 ,   0,   0]    ;%[0         0         0     ];
% Set them
h4_231(1).Color = black ; h4_232(1).Color = h4_231(1).Color; h4_233(1).Color = h4_231(1).Color;
h4_231(2).Color = yellow; h4_232(2).Color = h4_231(2).Color; h4_233(2).Color = h4_231(2).Color;
h4_231(3).Color = red   ; h4_232(3).Color = h4_231(3).Color; h4_233(3).Color = h4_231(3).Color;
h4_231(4).Color = blue  ; h4_232(4).Color = h4_231(4).Color; h4_233(4).Color = h4_231(4).Color;
% And set linestyles
h4_231(1).LineStyle = '-.' ; h4_232(1).LineStyle = h4_231(1).LineStyle; h4_233(1).LineStyle  = h4_231(1).LineStyle;
h4_231(2).LineStyle = ':' ; h4_232(2).LineStyle  = h4_231(2).LineStyle; h4_233(2).LineStyle  = h4_231(2).LineStyle;
h4_231(3).LineStyle = '--'; h4_232(3).LineStyle  = h4_231(3).LineStyle; h4_233(3).LineStyle  = h4_231(3).LineStyle;
h4_231(4).LineStyle = '-' ; h4_232(4).LineStyle  = h4_231(4).LineStyle; h4_233(4).LineStyle  = h4_231(4).LineStyle;

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-deps' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
savefig(        ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h', '.fig'])

%% Total Flow Field (u_tilde+u_bar|)
figure(825)
subplot(222)
u_bar_mesh = SF.u_bar_over_mesh(state.hk, state.rt, state.ue, state.ue / nu_inf);
surf(z_mesh/S, y_mesh/delta0, (state.u_tilde + u_bar_mesh )/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Predicted U/U_e   at x = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on
% Total Flow Field (Experimental)
subplot(221)
[z_exp_mesh, y_exp_mesh, u_exp_mesh] = ED.interp_u_mesh(state.x/h_VG);
surf(z_exp_mesh/S, y_exp_mesh/delta0, u_exp_mesh/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Experimental U/U_e   at X/h = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_fields_', num2str(state.x/h_VG), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.png'], '-r600');
savefig(        ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.fig'])
colormap(bone(16));
print('-deps' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.png'], '-r600');

%% And plot velocity profiles again for x=50h
% Set state for x=50h
new_state = state_50h;
state     = new_state;
% Filter and Compute spatial rate of change of u_tilde (reconstruction)
[state.u_tilde, u_bar, du_tilde_dx, dT_dx, mag_mesh] = ...
                    FM.u_tilde_rate_of_change( state, nu_inf);
% Open figure
figure(450)
% Centerline
subplot(231)
j_center = CM.j_center;
h4_231   = plot( new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_center )/ new_state.ue + new_state.u_tilde(:, j_center )/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, 0) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Central Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Mid Line
subplot(233)
j_quarter = round(j_center/2);
h4_233    = plot(new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_quarter)/ new_state.ue + new_state.u_tilde(:, j_quarter)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, -1) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
% title('Side Symmetry Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(-1)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% legend('Mixed Flow' , 'Pure Shear Flow', 'Total Predicted Flow', 'Experimental Result')
legend(h4_233([4 3 2 1]), 'Experiment' , 'Reconstruction', 'Shear flow', 'Interaction flow')
% Intermediate Line
subplot(232)
j_eighth = round((j_center+j_quarter)/2);
h4_232   = plot(new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    u_bar(:, j_eighth)/ new_state.ue + new_state.u_tilde(:, j_eighth)/ new_state.ue, CM.y_range/delta0, ...
    ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, - CM.z_range(j_eighth) ./ CM.S) / u_inf, CM.y_range/delta0, 'k');
grid on; axis([-0.3 1.1 0 1.5]); % grid on; axis([-0.5 1.5 0 2]);
% title('Intermediate Line');
title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0.5)])
xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
% Update colors
blue   = [0 , 147, 208]/255;%[0         0.4470    0.7410];
red    = [0 ,   0,   0]    ;%[0.8500    0.3250    0.0980];
yellow = [0 ,   0,   0]    ;%[0.9290    0.6940    0.1250];
black  = [0 ,   0,   0]    ;%[0         0         0     ];
% Set them
h4_231(1).Color = black ; h4_232(1).Color = h4_231(1).Color; h4_233(1).Color = h4_231(1).Color;
h4_231(2).Color = yellow; h4_232(2).Color = h4_231(2).Color; h4_233(2).Color = h4_231(2).Color;
h4_231(3).Color = red   ; h4_232(3).Color = h4_231(3).Color; h4_233(3).Color = h4_231(3).Color;
h4_231(4).Color = blue  ; h4_232(4).Color = h4_231(4).Color; h4_233(4).Color = h4_231(4).Color;
% And set linestyles
h4_231(1).LineStyle = '-.' ; h4_232(1).LineStyle = h4_231(1).LineStyle; h4_233(1).LineStyle  = h4_231(1).LineStyle;
h4_231(2).LineStyle = ':' ; h4_232(2).LineStyle  = h4_231(2).LineStyle; h4_233(2).LineStyle  = h4_231(2).LineStyle;
h4_231(3).LineStyle = '--'; h4_232(3).LineStyle  = h4_231(3).LineStyle; h4_233(3).LineStyle  = h4_231(3).LineStyle;
h4_231(4).LineStyle = '-' ; h4_232(4).LineStyle  = h4_231(4).LineStyle; h4_233(4).LineStyle  = h4_231(4).LineStyle;

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-deps' , ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
savefig(        ['fig/final_velocity_profiles_', num2str(state.x/h_VG)  , 'h', '.fig'])

%% Total Flow Field (u_tilde+u_bar|)
figure(850)
subplot(222)
u_bar_mesh = SF.u_bar_over_mesh(state.hk, state.rt, state.ue, state.ue / nu_inf);
surf(z_mesh/S, y_mesh/delta0, (state.u_tilde + u_bar_mesh )/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Predicted U/U_e   at x = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on
% Total Flow Field (Experimental)
subplot(221)
[z_exp_mesh, y_exp_mesh, u_exp_mesh] = ED.interp_u_mesh(state.x/h_VG);
surf(z_exp_mesh/S, y_exp_mesh/delta0, u_exp_mesh/ state.ue)
xlabel('Z/S'); ylabel('Y/\delta_{ref}'); colorbar;
title(['Experimental U/U_e   at X/h = ' num2str(round(state.x / h_VG*100)/100)]);
view(2); shading flat; caxis([0 1]);
axis([min(CM.z_range)/S max(CM.z_range)/S 0 2]);
box on

% Now print
set(gcf , 'PaperType', 'A5');
orient landscape;
print('-dpdf' , ['fig/final_velocity_fields_', num2str(state.x/h_VG), 'h'   , '.pdf']);
print('-depsc', ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.png'], '-r600');
savefig(        ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h'   , '.fig'])
colormap(bone(16));
print('-deps' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.eps']);
print('-dpng' , ['fig/final_velocity_fields_', num2str(state.x/h_VG)  , 'h_BW', '.png'], '-r600');