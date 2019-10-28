% Clean WorkSpace
clear all; close all; clc;

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

% Experimental Conditions
exp_cond.h     =  5.0e-3;                 % [m   ] height of vane
exp_cond.d     = 12.5e-3;                 % [m   ] distance between pair trailing edges
exp_cond.l     = 12.5e-3;                 % [m   ] vane chord lenght
exp_cond.D     = 30.0e-3;                 % [m   ] pair separation
exp_cond.AoA   = 18;                      % [m   ] pair separation
exp_cond.u_inf = 15.0   ;                 % [m/s ] edge velocity
% Define Crossflow scaling variables
S              = exp_cond.D / 2;          % [m   ] half-width of vortex system cell
h              = exp_cond.h;              % [m   ] height of vane
delta_ref      = 0.0150;                  % [m   ] ref boundary layer height (for gauge, need not be accurate!)
x              = exp_cond.h*d_list;       % [m   ] distance behing vane trailing edge

% % Read Data
% Read data
ED             = BLT_experiment_data_reader(fname_root, d_list, h);

% % Filter Data
% Set filters
ED.set_filters(d_list_flag, z_min_list, z_max_list, y_min_list, y_max_list);
% Filter data
ED.filter_data();

% % Set Scalling 
% Set Crossflow scaling variables
ED.set_crossflow_scaling_variables(S, delta_ref);

% % Now use the data!
% Try an interpolation
[u, v, w]      = ED.interp_velocities(24.999999, [0.1 0.2 0.3 0.4], 0);



target_y = linspace(0,3);
z_center =  0;
z_side   = -1;
u_center = ED.interp_velocities(5, target_y , z_center); 
u_side   = ED.interp_velocities(5, target_y , z_side  ); 
plot(u_center, target_y, u_side, target_y); grid on; hold on;
xlabel('Streamwise Speed'); ylabel('y/\delta - normal coordinate');






