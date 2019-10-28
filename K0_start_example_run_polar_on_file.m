%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, Ricardo Pereira, 2010-2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Run a polar on an airfoil file
%           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean environment
close all; clear all; clc;

%% Inputs
% % Filename
case_overview.filename       = 'naca0012.air';   % Name of airfoil file (can also have folders, absolute or from cd (./))
% % Polar conditions
case_overview.polar_range    = [-5 20 0.4 0];   % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
case_overview.Re             = 3000000   ;       % Reynolds Number (NaN means inviscid)
case_overview.xtr_top        =   0.99;           % set top    transition position for free
case_overview.xtr_bot        =   0.99;           % set bottom transition position for free
case_overview.N_crit         =   9   ;           % critical amplification factor for Tollmien-Schlichting waves

%% Initialize environment
% Add access to necessary paths
fs = filesep(); addpath([cd fs 'src']);

% Create and set system context
SC = system_context;
SC.N_cores = 1;
SC.set_context;

% Make simulation protocol object
SP1 = simulation_protocol('free_transition' , SC); SP1.target_application='RBINKOLIVEIRA_V2';
SP1.operation = 'alfa_polar_ref_start';                                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
SP1.operation_parameters = case_overview.polar_range;                           % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
SP1.Re                   = case_overview.Re;                                    % Reynolds Number (NaN means inviscid)
SP1.xtr_top              = case_overview.xtr_top;                          % set forced transition on top
SP1.xtr_bot              = case_overview.xtr_bot;                          % set forced transition on bottom
SP1.N_crit               = case_overview.N_crit;                                % critical amplification factor for Tollmien-Schlichting waves
SP1.save_bl_parameters   = 'true';

% Make simulation worker object
SW_1 = simulation_worker('simul_worker_1', SP1, [] , [],  SC);                      % Mariline
SW_1.app_name            = 'RBINKOLIVEIRA_V2';                                  % For xfoil write SW_1.app_name = 'xfoil'
SW_1.fig_active          = 0;                                                   % Activate plotting for clean case
SW_1.parallelized        = 1;


%% Now run polar
% Simulation worker runs everything and returns results in aerodynamic
% polar object
ap = SW_1.run_polar_on_airfoil_file(case_overview.filename);

%% Plot some results
% Using the plot method of the aerodynamic polar object
ap.plot();
% And on our own to plot boundary layer values
figure(2)
xun_range = linspace(0,1); aoa_plot = 5;
subplot(411)
plot(xun_range, ap.bldata.cpx_fun_top(aoa_plot, xun_range)); hold on;
plot(xun_range, ap.bldata.cpx_fun_bot(aoa_plot, xun_range)); grid on;
xlabel('x/c'); ylabel('C_p'); title('Pressure coefficient');
subplot(412)
plot(xun_range, ap.bldata.cfx_fun_top(aoa_plot, xun_range)); hold on;
plot(xun_range, ap.bldata.cfx_fun_bot(aoa_plot, xun_range));
xlabel('x/c'); ylabel('C_f'); title('Skin friction coefficient');
subplot(413)
plot(xun_range, ap.bldata.dst_fun_top(aoa_plot, xun_range)); hold on;
plot(xun_range, ap.bldata.dst_fun_bot(aoa_plot, xun_range));
xlabel('x/c'); ylabel('\delta^*/c'); title('Displacement thickness');
subplot(414)
plot(xun_range, ap.bldata.tet_fun_top(aoa_plot, xun_range)); hold on;
plot(xun_range, ap.bldata.tet_fun_bot(aoa_plot, xun_range));
xlabel('x/c'); ylabel('\theta/c'); title('Momentum thickness');


%% Accessing results in ordered arrays
% Monotonic array with angles of attack for which convergence was achieved
alpha_range = ap.alpha_range;
% Ordered Cl data
cl_range = ap.cl_alpha(alpha_range);
% Ordered Cd data
cd_range = ap.cd_alpha(alpha_range);
% Ordered Cm data
cm_range = ap.cm_alpha(alpha_range);
% Transition position on suction side
s_xtr_range = field_alpha(ap, 's_xtr', alpha_range);
% Transition position on pressure side
p_xtr_range = field_alpha(ap, 'p_xtr', alpha_range);
% Pressure drag
cdp_range = field_alpha(ap, 'cdp', alpha_range);




