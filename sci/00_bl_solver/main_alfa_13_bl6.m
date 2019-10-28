%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Integral Boundary Layer Integrator (ODE)
%           Plasma Development Tool
%
%       August 2014, GNU-GPLv3 or later
%       Gael de Oliveira, Ricardo Pereira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       This script:
%            1. Loads experimental data using the experimental_data_loader2
%               class 
%            2. Generates a case for the boundary_layer5_plasma ODE system
%               class using the plasma_descriptor class 
%            3. Integrates the BL equations with a very small fixed step
%               using the ode45 algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Startup

% Clear Stuff
clear all; clc; close all; 

% Add some useful paths
addpath 0_closure_relations
addpath 1_definitions
addpath 2_helper_functions

%% Define Inputs

% Choose experiment to reproduce
alpha_str='7';  % angle of attack considered
act_str='TE';  % Trailing or Leading Edge actuator
U_str='20';    % Wind Speed

% Ctau tuning parameter
ctau_ineq_factor = 1; 

% Starting point 
% (1 by default, greater to make sure transition is behind)
i_start = 12;


% Set Scales (semi-manual)
L           = 1;                                % [m]       Lenght Scale
rho_inf     = 1.225;                            % [kg/m3]   Free Stream Density
nu_inf      = 1.461e-5;                         % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
msq         = 0;                                % [--]      Unperturbed Free Stream Mach Number
U_inf       = str2double(U_str);                % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)
u_inf_over_nu_inf = U_inf / nu_inf;             % [1/m]


% Define CEI closure
dataset_file = 'cei_closure_data.mat';          % CEI Dataset filename

% Integration algorithm parameters
ode_solver  = @ode45;                          % ode45 or ode15s (stiff) are reasonable choices
fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps
N_steps     = 1000;                             % Set number of steps (only applicable for fixed_step mode


%% Load Experimental Data

%Generate filename strings
experiment_file_on=strcat('BL_',act_str,'_a',alpha_str,'_U',U_str,'_on.mat');
experiment_file_off=strcat('BL_',act_str,'_a',alpha_str,'_U',U_str,'_off.mat');

% Create EDL objects 
EDL_on = experimental_data_loader2();
EDL_off  = experimental_data_loader2();

% Set Scales and Case Parameters for Experiment
EDL_on.L           = L      ;  EDL_off.L           = L      ;
EDL_on.U_inf       = U_inf  ;  EDL_off.U_inf       = U_inf  ;
EDL_on.rho_inf     = rho_inf;  EDL_off.rho_inf     = rho_inf;
EDL_on.nu_inf      = nu_inf ;  EDL_off.nu_inf      = nu_inf ;


% Load data from experiment
EDL_on.load_experiment_data(experiment_file_on);
EDL_off.load_experiment_data(experiment_file_off);


%% Create Plasma Descriptor Objects

% Allocate PD_cell cell array
PD_cell = cell(1,EDL_on.N_p);

% Create a Plasma Descriptor for each Actuator
for n_p = 1:EDL_on.N_p
    % Create Object
    PD_cell{n_p} = plasma_descriptor();
    % Set Working Scales
    PD_cell{n_p}.set_adimensionalization_parameters(EDL_on.L, EDL_on.U_inf, EDL_on.rho_inf);
    % Set Inputs
    PD_cell{n_p}.set_plasma_dimensional_inputs( ...
        EDL_on.CFT_p_vector(n_p), ...
        EDL_on.L_p_vector(n_p), ...
        EDL_on.T_p_vector(n_p), ...
        EDL_on.X0_p_vector(n_p));    
end

%% Energy Interaction Coefficient (CEI) Closure

% Create manager object and load dataset (automatically)
CEI = cei_closure_manager(dataset_file);


%% Create Boundary Layer Object

% Create Boundary Layer Object
BL_on  = boundary_layer6_plasma();
BL_off = boundary_layer6();

% Set Hooks between BL object and and PD and CEI objects
BL_on.PD_cell   = PD_cell; 
BL_on.CEI  = CEI;

% Set nue 
BL_on.u_inf_over_nu_inf  = EDL_on.u_inf_over_nu_inf;
BL_off.u_inf_over_nu_inf = EDL_off.u_inf_over_nu_inf;

% Set forcing terms (order matters! must be done after u_inf_over_nu_inf was set)
msq_vector_on  = zeros(size(EDL_on.ue_vector));
BL_on.set_forcing_terms(EDL_on.xi_vector, EDL_on.ue_vector, msq_vector_on);

msq_vector_off = zeros(size(EDL_off.ue_vector));
BL_off.set_forcing_terms(EDL_off.xi_vector, EDL_off.ue_vector, msq_vector_off);



%% Prepare Integration

% % Extract Initial Conditions (Plasma on)
x0_on = EDL_on.xi_vector(i_start);                                  % Start of Integration 
t0_on = EDL_on.theta_vector(i_start);                               % Initial Momentum Thickness (adimensional=)
h0_on = EDL_on.dstr_vector(i_start) / EDL_on.theta_vector(i_start);          % Initial Shape Factor
% Set and guess ctau0 from equilibrium value
BL_on.set_equilibrium_initial_conditions(t0_on, h0_on);          % Initial Shear Stress Coefficient determined as if BL was in equilibrium (which more less never true!!!)

% % Extract Initial Conditions (Plasma off)
x0_off = EDL_off.xi_vector(i_start);                                  % Start of Integration 
t0_off = EDL_off.theta_vector(i_start);                               % Initial Momentum Thickness (adimensional=)
h0_off = EDL_off.dstr_vector(i_start) / EDL_off.theta_vector(i_start);          % Initial Shape Factor
% Set and guess ctau0 from equilibrium value
BL_off.set_equilibrium_initial_conditions(t0_off, h0_off);          % Initial Shear Stress Coefficient determined as if BL was in equilibrium (which more less never true!!!)

% % % Apply tuning correction for initial ctau 
% (really bad stuff... MUST REMAIN TEMPORARY)
BL_on.ctau0  = BL_on.ctau0  .* ctau_ineq_factor;
BL_off.ctau0 = BL_off.ctau0 .* ctau_ineq_factor;

% Determine end of integration (get either max of experimental data or
% plasma field actuation region)
xend_on     = max(max(EDL_on.X0_p_vector  + EDL_on.L_p_vector ), max(EDL_on.X )) / EDL_on.L ;
xend_off    = max(max(EDL_off.X0_p_vector + EDL_off.L_p_vector), max(EDL_off.X)) / EDL_off.L;

xend        = max(xend_on , xend_off);


%% Integration 
% Integrate the ODE!
BL_off.integrate_ODE(x0_off, xend);
BL_on.integrate_ODE( x0_on , xend);

%% Plot Results
% Extract some results
x_sol_on       = BL_on.sol.x      ; x_sol_off      = BL_off.sol.x      ;    
t_sol_on       = BL_on.sol.y(1,:) ; t_sol_off      = BL_off.sol.y(1,:) ;
h_sol_on       = BL_on.sol.y(2,:) ; h_sol_off      = BL_off.sol.y(2,:) ;
ctau_sol_on    = BL_on.sol.y(3,:) ; ctau_sol_off   = BL_off.sol.y(3,:) ;

% Retrieve dstr as it is more stable for comparison with experimental results!
%       Small (relative)oscilations in delta* (in the experimental results) 
%       lead to large oscillations of the shape factor (because delta* is small)
dstr_sol_on    = h_sol_on  .* t_sol_on;
dstr_sol_off   = h_sol_off .* t_sol_off;

% Plot the Shape Factor Evolution
figure(1)
subplot(2,1,1)
plot(x_sol_on , dstr_sol_on  , EDL_on.xi_vector , EDL_on.dstr_vector , ...
     x_sol_off, dstr_sol_off , EDL_off.xi_vector, EDL_off.dstr_vector);
legend('Plasma ON  - Computed', 'Plasma ON  - Experimental' , ...
       'Plasma OFF - Computed', 'Plasma OFF - Experimental' );
xlabel('x') ; ylabel('delta* - Displacement Thickness (adim.)') ; grid on

% Plot the Momentum Thickness Evolution
subplot(2,1,2)
plot(x_sol_on , t_sol_on  , EDL_on.xi_vector , EDL_on.theta_vector , ...
     x_sol_off, t_sol_off , EDL_off.xi_vector, EDL_off.theta_vector);
xlabel('x') ; ylabel('theta - Momentum Thickness (adim.)') ; grid on

figure(2)
subplot(2,1,1)
plot(x_sol_on , h_sol_on  , EDL_on.xi_vector , EDL_on.dstr_vector  ./ EDL_on.theta_vector , ...
     x_sol_off, h_sol_off , EDL_off.xi_vector, EDL_off.dstr_vector ./ EDL_off.theta_vector);
 legend('Plasma ON  - Computed', 'Plasma ON  - Experimental' , ...
       'Plasma OFF - Computed', 'Plasma OFF - Experimental' );
xlabel('x') ; ylabel('H12 - Shape Factor') ; grid on

subplot(2,1,2)
plot(EDL_on.xi_vector , EDL_on.ue_vector, EDL_off.xi_vector , EDL_off.ue_vector)
xlabel('x') ; ylabel('ue - Edge Speed (adim.)') ; grid on


