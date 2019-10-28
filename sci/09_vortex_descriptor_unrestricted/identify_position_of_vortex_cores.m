
% % Header copied from:
%       sct_passive_scalar_convection_alfa4_3AF_picture_mod_alfa38ZPG

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
x0        =    0 * h_VG;              % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
x_end     =    50.25* h_VG;%0.25;%;0.5;        % End Time (0.25 = 50 * h_vg)
dx        =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
dx_info   =    0.25 * h_VG;           % Time between information events
information_period=round(dx_info/dx); % Steps between information events

% Define output parameters
user_prompt = false;
img_output = true;
mat_output = true;
baldwin_factor = 1;

% Set new step size
pod_dx                 = 1.25e-3;     % Up to 1000 times bigger than the dx of the normal stuff!
% Set pod information period
pod_information_period = round(dx_info/pod_dx);


% Set Outer Flow Distribution Parameters
xi_vector  = linspace(x0,x_end,100);            % Reinterpolated, number of points (100) not so relevant
ue_vector  = u_inf   * ones(size(xi_vector));   % Set flat pressure gradient with u_inf outer speed
msq_vector = msq_inf * ones(size(xi_vector));   % Set all flow as purely incompressible

%% Initialization: Create Vortex System
n_cells = 5;
S       = D_VG/2;

VD = vortex_descriptor(S, n_cells);
VD.induction_function = 3;              % Try Lamb Vortex

%% Initializaion: Make z-y Meshes

CM = crossflow_mesh(S, delta0);

z_mesh = CM.z_mesh;
y_mesh = CM.y_mesh;



%% Now load experimental results
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


%% Identify location of vortex cores with plots! (quiver plot works best!)
figure(1)
subplot(221)
n = 1; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat
title('|V+W|')
subplot(222)
n = 6; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat
subplot(223)
n = 7; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat
subplot(224)
n = 8; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat

figure(2)

subplot(221)
n = 1; 
[curlx, cav] = curl(ED.data_struct_cell{n}.w,ED.data_struct_cell{n}.v);
surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, curlx); view(2); shading flat;
caxis([-1, 1]*0.1)
title('\omega_x')
subplot(222)
n = 2; 
[curlx, cav] = curl(ED.data_struct_cell{n}.w,ED.data_struct_cell{n}.v);
surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, curlx); view(2); shading flat;
caxis([-1, 1]*0.1)
subplot(223)
n = 3; 
[curlx, cav] = curl(ED.data_struct_cell{n}.w,ED.data_struct_cell{n}.v);
surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, curlx); view(2); shading flat;
caxis([-1, 1]*0.1)
subplot(224)
n = 4; 
[curlx, cav] = curl(ED.data_struct_cell{n}.w,ED.data_struct_cell{n}.v);
surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, curlx); view(2); shading flat;
caxis([-1, 1]*0.1)


figure(3)
subplot(221)
n = 1; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.v); view(2); shading flat
title('V')
subplot(222)
n = 6; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.v); view(2); shading flat
subplot(223)
n = 7; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.v); view(2); shading flat
subplot(224)
n = 8; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.v); view(2); shading flat


figure(4)
subplot(221)
n = 1; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat
title('|V+W|')
subplot(222)
n = 6; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat
subplot(223)
n = 7; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat
subplot(224)
n = 8; surf(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, sqrt(ED.data_struct_cell{n}.v.^2 + ED.data_struct_cell{n}.w.^2)); view(2); shading flat


% % New approach
figure(5)
n = 1;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);


figure(6)
n = 2;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);

figure(7)
n = 3;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);

figure(8)
n = 4;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);

figure(9)
n = 5;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);

figure(10)
n = 6;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);

figure(25)
n = 7;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);

figure(50)
n = 8;
quiver(ED.data_struct_cell{n}.z, ED.data_struct_cell{n}.y, ED.data_struct_cell{n}.w, ED.data_struct_cell{n}.v);

%% And write them down!

% % Move on to write down what we just observed/measured manually
y_l     = [ 3.40 , 3.45 , 3.61 , 3.72, 3.88, 4.10] * 1e-3;
z_l     = [-6.45 ,-6.95 ,-7.28 ,-7.88,-8.30,-8.48] * 1e-3;
y_r     = [ 3.08 , 3.15 , 3.12 , 3.12, 3.20, 3.25] * 1e-3;
z_r     = [ 6.70 , 7.05 , 7.65 , 7.95, 8.35, 8.61] * 1e-3;
d_list  = [ 5    , 6    , 7    , 8   , 9   , 10  ]       ;

% % Now make a similar one for longer stuff
y_lB    = [ 3.40 , 3.45 , 3.61 , 3.72, 3.88, 4.10, 6.6 ] * 1e-3;
z_lB    = [-6.45 ,-6.95 ,-7.28 ,-7.88,-8.30,-8.48,-10.1] * 1e-3;
y_rB    = [ 3.08 , 3.15 , 3.12 , 3.12, 3.20, 3.25, 4.6 ] * 1e-3;
z_rB    = [ 6.70 , 7.05 , 7.65 , 7.95, 8.35, 8.61, 10.5] * 1e-3;
d_listB = [ 5    , 6    , 7    , 8   , 9   , 10  , 25  ];


%% Now plot, make symmetry average and error bars based on asymmetry!

% Now plot
figure(101)
plot(d_list, -z_l, d_list, z_r); grid on; legend('Left (from back)', 'Right (from back)')
figure(102)
plot(d_list,  y_l, d_list, y_r); grid on; legend('Left (from back)', 'Right (from back)')

% Now reprocess list of values into average
y_mean  = 0.5 * ( y_l + y_r);
z_mean  = 0.5 * (-z_l + z_r);

y_meanB  = 0.5 * ( y_lB + y_rB);
z_meanB  = 0.5 * (-z_lB + z_rB);

% And unbiased standard deviation (Bessel correction)
std_y   = sqrt( (( y_l - y_mean).^2 + (y_r - y_mean).^2 ) / (2 - 1) );
std_z   = sqrt( ((-z_l - z_mean).^2 + (z_r - z_mean).^2 ) / (2 - 1) );

std_yB  = sqrt( (( y_lB - y_meanB).^2 + (y_rB - y_meanB).^2 ) / (2 - 1) );
std_zB  = sqrt( ((-z_lB - z_meanB).^2 + (z_rB - z_meanB).^2 ) / (2 - 1) );

% Now replot all with error bars
figure(101)
hold on;
errorbar(d_list, z_mean, std_z);

figure(102)
hold on;
errorbar(d_list, y_mean, std_y);


% Now make offsets for core heights
y_offset_r = zeros(size(d_list));
y_offset_l = zeros(size(d_list));
% One by one
for n = 1:6
    y_offset_l(n) = interp1(ED.data_struct_cell{n}.z(end,:), ED.data_struct_cell{n}.y_offset_extrapolation, z_l(n));
    y_offset_r(n) = interp1(ED.data_struct_cell{n}.z(end,:), ED.data_struct_cell{n}.y_offset_extrapolation, z_r(n));
end
% And unbias y
y_lu = y_l + y_offset_l;
y_ru = y_r + y_offset_r;

% Now make offsets for core heights on longer one
y_offset_rB = zeros(size(d_listB));
y_offset_lB = zeros(size(d_listB));
% One by one
for n = 1:7
    y_offset_lB(n) = interp1(ED.data_struct_cell{n}.z(end,:), ED.data_struct_cell{n}.y_offset_extrapolation, z_lB(n));
    y_offset_rB(n) = interp1(ED.data_struct_cell{n}.z(end,:), ED.data_struct_cell{n}.y_offset_extrapolation, z_rB(n));
end
% And unbias y
y_luB = y_lB + y_offset_lB;
y_ruB = y_rB + y_offset_rB;



figure(102)
hold on;
plot(d_list, y_lu, 'k--')
plot(d_list, y_ru, 'k--')

% So, now redo averages
yu_mean  = 0.5 * ( y_lu + y_ru);
yu_meanB = 0.5 * ( y_luB + y_ruB);

% So, now redo standard deviation
std_yu  = sqrt( (( y_lu - yu_mean).^2 + (y_ru - yu_mean).^2 ) / (2 - 1) );
std_yuB  = sqrt( (( y_luB - yu_meanB).^2 + (y_ruB - yu_meanB).^2 ) / (2 - 1) );

figure(102)
hold on;
errorbar(d_list, yu_mean, std_yu);

%% Now load some numerical results
% Now get numerical results manually
y_num     = zeros(size(yu_mean));
z_num     = zeros(size(z_mean));
state_05h = load('20180727T013559/matfiles/full_state_asy_5h.mat' ); y_num(1) = state_05h.state.y_v; z_num(1) = state_05h.state.z_v;
state_06h = load('20180727T013559/matfiles/full_state_asy_6h.mat' ); y_num(2) = state_06h.state.y_v; z_num(2) = state_06h.state.z_v;
state_07h = load('20180727T013559/matfiles/full_state_asy_7h.mat' ); y_num(3) = state_07h.state.y_v; z_num(3) = state_07h.state.z_v;
state_08h = load('20180727T013559/matfiles/full_state_asy_8h.mat' ); y_num(4) = state_08h.state.y_v; z_num(4) = state_08h.state.z_v;
state_09h = load('20180727T013559/matfiles/full_state_asy_9h.mat' ); y_num(5) = state_09h.state.y_v; z_num(5) = state_09h.state.z_v;
state_10h = load('20180727T013559/matfiles/full_state_asy_10h.mat'); y_num(6) = state_10h.state.y_v; z_num(6) = state_10h.state.z_v;
% And (semi)automatically
y_num_full = zeros(1,49);
z_num_full = zeros(1,49);
for n_num=1:length(y_num_full)
   state_tmp = load(['20180727T013559/matfiles/full_state_asy_', num2str(n_num), 'h.mat' ]); 
   y_num_full(n_num) = state_tmp.state.y_v; 
   z_num_full(n_num) = state_tmp.state.z_v; 
end
state_tmp = load(['20180727T013559/matfiles/full_state_asy_', '-0.75', 'h.mat' ]); 
y_num_full = [state_tmp.state.y_v , y_num_full];
z_num_full = [state_tmp.state.z_v , z_num_full];

% And plot comparison
figure(301)
errorbar(d_list, z_mean, std_z); hold on; grid on;
plot(d_list, z_num)

figure(302)
errorbar(d_list, yu_mean, std_yu); hold on; grid on;
plot(d_list, y_num)


% Now make joint representation (not so nice)
figure(303)
plot( z_l  , y_lu, 'x-k'); hold on
plot( z_r  , y_ru, 'x-k');
plot(-z_num, y_num, 'x-r')
plot( z_num, y_num, 'x-r')

% And comparison plot
figure(304);

subplot(2,2,[1 2])
d_listB_subset = [1 6 7];
% Plot left side
h0 = plot(z_lB(d_listB_subset) / S, y_luB(d_listB_subset) / delta0, 'o-'); hold on;
h1 = plot(-z_num_full / S, y_num_full /delta0 , '.-'); 
plot(  1     * (d_VG / 2) / S,    1  * h_VG /delta0, 'k-x')
plot( [1, 1] * (d_VG / 2) / S, [0 1] * h_VG /delta0, 'k-')

% Plot right side
h2 = plot(z_rB(d_listB_subset) / S, y_ruB(d_listB_subset) / delta0, 'o-');
h3 = plot( z_num_full / S, y_num_full /delta0 , '.-'); grid on;
plot(- 1     * (d_VG / 2) / S,    1  * h_VG /delta0, 'k-x')
plot(-[1, 1] * (d_VG / 2) / S, [0 1] * h_VG /delta0, 'k-')

% Make symmetry line
plot([0 0], [0,0.8], 'k-.')
% Reset axis
axis([-0.8000    0.8000         0    0.8000])
% Agregate colors
h2.Color = h0.Color;
h3.Color = h1.Color;
% Decorate!
xlabel('z/S'); ylabel('y/\delta_{ref}');
legend('Experimental', 'Prediction (x = 0 to 50h)', 'Vane Trailing Edge', 'Location', 'North')
title('Crossflow movement of vortex filament core')


% Add labels to experimental data
nB_subset=1; text(z_rB(d_listB_subset(nB_subset)) / S, y_rB(d_listB_subset(nB_subset)) / delta0 - 0.02, '5h')
nB_subset=2; text(z_rB(d_listB_subset(nB_subset)) / S, y_rB(d_listB_subset(nB_subset)) / delta0 - 0.02, '10h')
nB_subset=3; text(z_rB(d_listB_subset(nB_subset)) / S, y_rB(d_listB_subset(nB_subset)) / delta0 - 0.02, '25h')

nB_subset=1; text(z_lB(d_listB_subset(nB_subset)) / S - 0.05, y_lB(d_listB_subset(nB_subset)) / delta0 - 0.02, '5h')
nB_subset=2; text(z_lB(d_listB_subset(nB_subset)) / S - 0.08, y_lB(d_listB_subset(nB_subset)) / delta0 - 0.02, '10h')
nB_subset=3; text(z_lB(d_listB_subset(nB_subset)) / S + 0.02, y_lB(d_listB_subset(nB_subset)) / delta0 + 0.05, '25h')

subplot(223)
errorbar(d_list, z_mean/S, 2 * mean(std_z) * ones(size(std_z)) /S ); hold on
plot(d_list, z_num / S, '.-'); grid on;
xlabel('x/h^{VG}'); ylabel('z/S'); xlabel('x/h');
legend('Experimental', 'Prediction', 'Location', 'SouthEast')
title('Spanwise movement of vortex core')
axis([5.0000   10.0000    0.4000    0.6000]);

subplot(224)
errorbar(d_list, y_mean/delta0, 2 * mean(std_y) * ones(size(std_y)) /delta0 ); hold on
plot(d_list, y_num / S, '.-'); grid on;
xlabel('x/h^{VG}'); ylabel('y/\delta_{ref}'); xlabel('x/h');
legend('Experimental', 'Prediction', 'Location', 'SouthEast')
title('Normal movement of vortex core')
axis([5.0000   10.0000    0.1500    0.3500]);

set(gcf , 'PaperType', 'A5');
orient landscape;
print ('-dpdf', ['vortex_core_positions_ZPG_extrapolated.pdf']) %#ok<NBRAK>
savefig vortex_core_positions_ZPG_extrapolated.fig


%% Now compare induction of vortex system along a vertical line

% Choose stance (1 = 5h)
nd_list = 2;
state   = state_06h.state;

% Make comparison line
y_line     = linspace(0,4);
z_line_exp =   z_l(nd_list ) * ones(size(y_line)) / S;
z_line_num = - state.z_v     * ones(size(y_line)) / S;

% Get experimental spanwise velocity over comparison lines
w_line = interp2(ED.data_struct_cell{nd_list}.z / S, ED.data_struct_cell{nd_list}.y / delta0, ED.data_struct_cell{nd_list}.w, z_line_exp , y_line);


% Now compute vortical field prediction
[w_mesh, v_mesh, mag_mesh] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v);
% Get numerical prediciotn for spanwise velocity over comparison line
w_line_num = interp2(CM.z_mesh / S, CM.y_mesh / delta0, w_mesh, z_line_num , y_line);


% % Now plot

% First plot raw
figure(5)
plot(y_line, w_line, y_line, w_line_num); hold on; grid on
% Then add smoothed and shifted lines (this is LEGIT: there is something
% wrong with the terminal w velocity)
plot(y_line, smooth(w_line), y_line, smooth(smooth(w_line)))
grid on

% Then smooth 
w_line_smoothed = smooth(smooth(smooth(w_line)));

% And offset (to correct for systematic error on far-field spanwise
% velocity)
w_offset = 0.11;
w_line_smoothed_offsetted = w_line_smoothed + w_offset;
w_line_offsetted          = w_line          + w_offset;


figure(6)
plot(y_line, w_line_smoothed_offsetted ); hold on;
plot(y_line,               w_line_num  ); grid on;
axis([0 2 -0.05 0.05])
xlabel('z/S'); ylabel('w/U_e');
title('Spanwise induction of vortex core at x = 6h')


%% New approach to core visualization! No smoothing needed!
% Choose stance (2 = 6h)
nd_list = 2;
state   = state_06h.state;
% Make comparison line
y_line     = linspace(0,4);
z_line_exp =   z_l(nd_list ) * ones(size(y_line)) / S;
z_line_num = - state.z_v     * ones(size(y_line)) / S;

% Get experimental spanwise velocity over comparison lines
w_line06 = interp2(ED.data_struct_cell{nd_list}.z / S, ED.data_struct_cell{nd_list}.y / delta0, ED.data_struct_cell{nd_list}.w, z_line_exp , y_line);
% Now compute vortical field prediction
[w_mesh, ~, ~] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v/2);
% Get numerical prediciotn for spanwise velocity over comparison line
w_line_num06 = interp2(CM.z_mesh / S, CM.y_mesh / delta0, w_mesh, z_line_num , y_line);


nd_list = 4;
state   = state_08h.state;
% Make comparison line
y_line     = linspace(0,4);
z_line_exp =   z_l(nd_list ) * ones(size(y_line)) / S;
z_line_num = - state.z_v     * ones(size(y_line)) / S;

% Get experimental spanwise velocity over comparison lines
w_line08 = interp2(ED.data_struct_cell{nd_list}.z / S, ED.data_struct_cell{nd_list}.y / delta0, ED.data_struct_cell{nd_list}.w, z_line_exp , y_line);
% Now compute vortical field prediction
[w_mesh, v_mesh, mag_mesh] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v/2);
% Get numerical prediciotn for spanwise velocity over comparison line
w_line_num08 = interp2(CM.z_mesh / S, CM.y_mesh / delta0, w_mesh, z_line_num , y_line);


nd_list = 6;
state   = state_10h.state;
% Make comparison line
y_line     = linspace(0,4);
z_line_exp =   z_l(nd_list ) * ones(size(y_line)) / S;
z_line_num = - state.z_v     * ones(size(y_line)) / S;

% Get experimental spanwise velocity over comparison lines
w_line10 = interp2(ED.data_struct_cell{nd_list}.z / S, ED.data_struct_cell{nd_list}.y / delta0, ED.data_struct_cell{nd_list}.w, z_line_exp , y_line);
% Now compute vortical field prediction
[w_mesh, v_mesh, mag_mesh] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v/2);
% Get numerical prediciotn for spanwise velocity over comparison line
w_line_num10 = interp2(CM.z_mesh / S, CM.y_mesh / delta0, w_mesh, z_line_num , y_line);

figure(106)
h_06 = plot(y_line, w_line06/state.ue); hold on
h_08 = plot(y_line, w_line08/state.ue);
h_10 = plot(y_line, w_line10/state.ue);

h_06 = plot(y_line, w_line_num06/state.ue);
h_08 = plot(y_line, w_line_num08/state.ue);
h_10 = plot(y_line, w_line_num10/state.ue);

legend('6', '8', '10', '6n', '8n', '10n')

figure(206)

h_06 = plot(y_line, w_line06/state.ue); hold on
h_10 = plot(y_line, w_line10/state.ue);

h_06 = plot(y_line, w_line_num06/state.ue);
h_10 = plot(y_line, w_line_num10/state.ue);

legend('6', '10', '6n', '10n')


%% Yet a new approach to core visualization! No fiddling needed!
% Choose stance (2 = 6h)
nd_list = 2;
state   = state_06h.state;
% Make comparison line
y_line     = linspace(0,4);
z_line_exp =   z_l(nd_list ) * ones(size(y_line)) / S;
z_line_num = z_line_exp; %- state.z_v     * ones(size(y_line)) / S;

% Get experimental spanwise velocity over comparison lines
w_line06 = interp2(ED.data_struct_cell{nd_list}.z / S, ED.data_struct_cell{nd_list}.y / delta0, ED.data_struct_cell{nd_list}.w, z_line_exp , y_line);
% Now compute vortical field prediction
[w_mesh, ~, ~] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v/2);
% Get numerical prediciotn for spanwise velocity over comparison line
w_line_num06 = interp2(CM.z_mesh / S, CM.y_mesh / delta0, w_mesh, z_line_num , y_line);


nd_list = 6;
state   = state_10h.state;
% Make comparison line
y_line     = linspace(0,4);
z_line_exp =   z_l(nd_list ) * ones(size(y_line)) / S;
z_line_num =   z_line_num; % - state.z_v     * ones(size(y_line)) / S;

% Get experimental spanwise velocity over comparison lines
w_line10 = interp2(ED.data_struct_cell{nd_list}.z / S, ED.data_struct_cell{nd_list}.y / delta0, ED.data_struct_cell{nd_list}.w, z_line_exp , y_line);
% Now compute vortical field prediction
[w_mesh, v_mesh, mag_mesh] = VD.induction_on_mesh(CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v/2);
% Get numerical prediciotn for spanwise velocity over comparison line
w_line_num10 = interp2(CM.z_mesh / S, CM.y_mesh / delta0, w_mesh, z_line_num , y_line);


figure(306)

h_06 = plot(y_line, w_line06/state.ue); hold on
h_10 = plot(y_line, w_line10/state.ue);

h_06 = plot(y_line, w_line_num06/state.ue);
h_10 = plot(y_line, w_line_num10/state.ue);

legend('6', '10', '6n', '10n')
