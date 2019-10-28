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
clear all; close all; clc; 

% Add some useful paths
addpath closure_relations
addpath definitions

%% Define Inputs

% Ctau tuning parameter
ctau_ineq_factor = 1; 

% Set Scales (semi-manual)
L           = 1;                                % [m]       Lenght Scale
rho_inf     = 1.225;                            % [kg/m3]   Free Stream Density
nu_inf      = 1.461e-5;                         % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
msq_inf     = 0;                                % [--]      Unperturbed Free Stream Mach Number
u_inf       = 15;                               % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)
u_inf_over_nu_inf = u_inf / nu_inf;             % [1/m]

% Integration algorithm parameters
ode_solver  = @ode45;                           % ode45 or ode15s (stiff) are reasonable choices
fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps
N_steps     = 1000;                             % Set number of steps (only applicable for fixed_step mode

% External forcing parameters
xi_vector  = linspace(0,20,100);                % 
ue_vector  = u_inf   * ones(size(xi_vector));   % Set flat pressure gradient with u_inf outer speed
msq_vector = msq_inf * ones(size(xi_vector));   % Set all flow as purely incompressible

% Initial Conditions
x0   = xi_vector(1);            % Start of Integration 
rt0  = 500;                    % Initial Momentum Thickness Reynolds
t0   = rt0/u_inf_over_nu_inf;   % Initial Momentum Thickness (adimensional=)
h0   = 1.4;                     % Initial Shape Factor

%% Create Boundary Layer Object

% Create Boundary Layer Object
BL_off = boundary_layer6();

% Set nue 
BL_off.u_inf_over_nu_inf = u_inf_over_nu_inf;

% Set forcing terms (order matters! must be done after u_inf_over_nu_inf was set)
BL_off.set_forcing_terms(xi_vector, ue_vector, msq_vector);


%% Prepare Integration

% % Extract Initial Conditions


% Set and guess ctau0 from equilibrium value
BL_off.set_equilibrium_initial_conditions(t0, h0);          % Initial Shear Stress Coefficient determined as if BL was in equilibrium (which more less never true!!!)

% % % Apply tuning correction for initial ctau 
% (really bad stuff... MUST REMAIN TEMPORARY, does nothing for ctau_ineq_factor=1)
BL_off.ctau0 = BL_off.ctau0 .* ctau_ineq_factor;

% Determine end of integration
xend        = xi_vector(end);


%% Integration 
% Integrate the ODE!
BL_off.integrate_ODE(x0, xend);

%% Plot Results
% Extract some results
x_sol_off      = BL_off.sol.x      ; 
t_sol_off      = BL_off.sol.y(1,:) ; 
h_sol_off      = BL_off.sol.y(2,:) ; 
ctau_sol_off   = BL_off.sol.y(3,:) ; 

% Retrieve dstr as it is more stable for comparison with experimental results!
%       Small (relative)oscilations in delta* (in the experimental results) 
%       lead to large oscillations of the shape factor (because delta* is small)
dstr_sol_off   = h_sol_off .* t_sol_off;

% Plot the Shape Factor Evolution
figure(1)
subplot(2,1,1)
plot(x_sol_off, dstr_sol_off);
grid on;
xlabel('x') ; ylabel('delta* - Displacement Thickness (adim.)') ; 

% Plot the Momentum Thickness Evolution
subplot(2,1,2)
plot(x_sol_off, t_sol_off );
xlabel('x') ; ylabel('theta - Momentum Thickness (adim.)') ; grid on

% Save to file
savefig(['figures/BL1_',datestr(now, 30),'.fig']);


figure(2)
subplot(2,1,1)
plot(x_sol_off, h_sol_off);
xlabel('x') ; ylabel('H12 - Shape Factor') ; grid on

subplot(2,1,2)
plot(xi_vector , ue_vector)
xlabel('x') ; ylabel('ue - Edge Speed (adim.)') ; grid on

% Save to file
savefig(['figures/BL2_',datestr(now, 30),'.fig']);

%% Write Comparison Strings

disp(['Case              |      dstr       |        theta        |'])
disp(['Clean             |    ' num2str(dstr_sol_off(end)) , '   |      ' , num2str(t_sol_off(end)) , '     |']);







