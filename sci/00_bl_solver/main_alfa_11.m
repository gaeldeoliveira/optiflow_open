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
experiment_file = 'BL_LE_a13_U20_on.mat';

% Set Scales 
L           = 1;                                % [m]       Lenght Scale
U_inf       = 20;                               % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)
rho_inf     = 1.225;                            % [kg/m3]   Free Stream Density
nu_inf      = 1.461e-5;                         % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
u_inf_over_nu_inf = U_inf / nu_inf;             % [1/m]
msq         = 0;                                % [--]      Unperturbed Free Stream Mach Number

% Define CEI closure
dataset_file = 'cei_closure_data.mat';          % CEI Dataset filename

% Integration algorithm parameters
ode_solver  = @ode15s;                          % ode45 or ode15s (stiff) are reasonable choices
fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps
N_steps     = 1000;                             % Set number of steps (only applicable for fixed_step mode



%% Load Experimental Data

% Create EDL object 
EDL = experimental_data_loader2();

% Set Scales and Case Parameters for Experiment
EDL.L           = L;
EDL.U_inf       = U_inf;
EDL.rho_inf     = rho_inf;
EDL.nu_inf      = nu_inf;

% Load data from experiment
EDL.load_experiment_data(experiment_file);


%% Create Plasma Descriptor Objects

% Allocate PD_cell cell array
PD_cell = cell(1,EDL.N_p);

% Create a Plasma Descriptor for each Actuator
for n_p = 1:EDL.N_p
    % Create Object
    PD_cell{n_p} = plasma_descriptor();
    % Set Working Scales
    PD_cell{n_p}.set_adimensionalization_parameters(EDL.L, EDL.U_inf, EDL.rho_inf);
    % Set Inputs
    PD_cell{n_p}.set_plasma_dimensional_inputs( ...
        EDL.CFT_p_vector(n_p), ...
        EDL.L_p_vector(n_p), ...
        EDL.T_p_vector(n_p), ...
        EDL.X0_p_vector(n_p));    
end

%% Energy Interaction Coefficient (CEI) Closure

% Create manager object and load dataset (automatically)
CEI = cei_closure_manager(dataset_file);


%% Create Boundary Layer Object

% Create Boundary Layer Object
BL = boundary_layer5_plasma();

% Set Hooks between BL object and and PD and CEI objects
BL.PD_cell   = PD_cell; 
BL.CEI  = CEI;

% Set nue 
BL.u_inf_over_nu_inf = EDL.u_inf_over_nu_inf;

% Set forcing terms (order matters! must be done after u_inf_over_nu_inf was set)
msq_vector = zeros(size(EDL.ue_vector));
BL.set_forcing_terms(EDL.xi_vector, EDL.ue_vector, msq_vector);


%% Prepare Integration

% Determine start of integration 
x0    = EDL.xi_vector(1); 

% Determine end of integration (get either max of experimental data or
% plasma field actuation region)
Xend = max(max(EDL.X0_p_vector + EDL.L_p_vector), max(EDL.X));
xend  = Xend / EDL.L;

% % Extract Initial Conditions
x = EDL.xi_vector(1);
t = EDL.theta_vector(1);
h = EDL.dstr_vector(1) / EDL.theta_vector(1);

% % Guess corresponding ctau (assume equilibrium ctau)
% From equilibrium ctau (even though even a flat plate is not in perfect equilibrium!)
hk      = hkin( h, msq);                    
rt      = re_theta( t, BL.u_inf_over_nu_inf ); % Looks ok
cf      = cft( hk, rt, msq);                % Looks ok
hs      = hst( hk, rt, msq);
us      = usg(h, hk, hs);
ctau    = ctauzero( h , hk , hs , us);      % Looks ok

% % Set initial conditions 
y0(1) = t;
y0(2) = h;
y0(3) = ctau;

% Define handle function reference to wrapper
odefun = @BL.ode_fun_wrapper;

%% Integration 

% Integrate
if fixed_step  == false
    % Variable Step Size
    sol = ode_solver(odefun,[x0 xend],y0);
else
    % Fixed Step Size
    sol = ode_solver(odefun,linspace(x0, xend, N_steps),y0);
end

%% Plot Results
% Extract some results
t_sol       = sol.y(1,:);
h_sol       = sol.y(2,:);
ctau_sol    = sol.y(3,:);

% Plot the Shape Factor Evolution
figure(1)
plot(sol.x, h_sol .* t_sol , EDL.xi_vector , EDL.dstr_vector)
legend('Computed Solution', 'Experimental Result')
xlabel('x') ; ylabel('delta* - Displacement Thickness') ; grid on

% Plot the Momentum Thickness Evolution
figure(2)
plot(sol.x, t_sol , EDL.xi_vector , EDL.theta_vector)
legend('Computed Solution', 'Experimental Result')
xlabel('x') ; ylabel('theta - Momentum Thickness') ; grid on

