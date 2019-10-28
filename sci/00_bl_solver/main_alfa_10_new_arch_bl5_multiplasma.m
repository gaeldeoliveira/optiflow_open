%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Integral Boundary Layer Integrator (ODE)
%           Plasma Development Tool
%
%       August 2014, GNU-GPLv3 or later
%       Gael de Oliveira, Ricardo Pereira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Todo List
%   Insert closure relations from Rfoilsuc code
%   Make velocity profile generator
%   Make agregator module as object, for plotting and closure set definition (using handle functions)
%   Make integrator module as object

% Clear Stuff
clear all; clc; % close all; 


% Add some useful paths
addpath 0_closure_relations
addpath 1_definitions
addpath 2_helper_functions

%% Inputs
% Set Scales 
L                   = 1;                        % [m]   Lenght Scale (Longitudinal)
U_inf               = 20;                       % [m/s] Unperturbed Flow Speed
rho_inf             = 1.225;                    % [kg/m3]Unperturbed Flow Density
u_inf_over_nu_inf   = 1e6;                      % [1/m] 

% Define plasma Parameters
F_T                 = -0.2;                     % [N/m] Total Actuator Force per Unit Span
L_p                 = 0.0975;                   % [m]   Plasma Force Field Lenght
T_p                 = 0.003;                    % [m]   Plasma Force Field Thickness
X0_p                = 0.2;                      % [m]   Plasma Force Field Start

% Define Boundary Layer Initial Conditions
X0                  = 0.0;                      % [m]   Starting Point (dimensional, will be adimensionalized to longitudinal lenght scale, no translation)
delta_2             = 0.008;                    % [m]   Momentum Thickness (dimensional, will be adimensionalized to longitudinal lenght scale)
H12                 = 1.36;                     % [--]  Shape Factor 

% Define Forcing Terms
xi_vector   = linspace(0, 0.5, 100);            % xi  = X  / L          Adimensional Stance       
ue_vector   = ones(size(xi_vector));            % ue  = Ue / Uinf       Adimensional Edge Speed
msq_vector  = zeros(size(xi_vector));           % msq = Ue / c          Edge Mach Number

% Define CEI closure
dataset_file = 'cei_closure_data.mat';          % CEI Dataset filename

% Integration algorithm parameters
ode_solver  = @ode15s;                          % ode45 or ode15s (stiff) are reasonable choices
fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps
N_steps     = 1000;                             % Set number of steps (only applicable for fixed_step mode


%% Allocation and Instanciation
% % % Energy Interaction Coefficient (CEI) Object
% Create object and load dataset (automatically)
CEI = cei_closure_manager(dataset_file);

% % % Plasma Descriptor Object
% Create Object
PD = plasma_descriptor();
% Set Working Scales
set_adimensionalization_parameters(PD, L, U_inf, rho_inf);
% Set Inputs
PD.set_plasma_dimensional_inputs(F_T, L_p, T_p, X0_p)

% % % Boundary Layer Object
% Create Boundary Layer Object
BL = boundary_layer5_plasma();
% Set Hooks between BL object and and PD and CEI objects
BL.PD   = PD; BL.CEI  = CEI;
% Set nue 
BL.u_inf_over_nu_inf = u_inf_over_nu_inf;
% Set forcing terms (order matters! must be done after u_inf_over_nu_inf)
BL.set_forcing_terms(xi_vector, ue_vector, msq_vector);


%% Build Consistent Initial Conditions

% % Adimensionalize Initial Conditions (ctau will be guessed later)
x       = 0.0   / L;                            % [--]  Starting Stance (adimensionalized to longitudinal lenght scale)
t       = 0.003 / L;                            % [--]  Momentum Thickness (adimensionalized to longitudinal lenght scale)
h       = H12;                                  % [--]  Shape Factor
msq     = 0;                                    % [--]  Mach Number


% % Guess corresponding ctau 
% From equilibrium ctau (even though flat plate is not in perfect equilibrium!)
hk      = hkin( h, msq);                    
rt      = re_theta( t, u_inf_over_nu_inf );       % Looks ok
cf      = cft( hk, rt, msq);                % Looks ok
hs      = hst( hk, rt, msq);
us      = usg(h, hk, hs);
ctau    = ctauzero( h , hk , hs , us);      % Looks ok
cdi     = cdissipation( us, cf, ctau);      % Looks ok


%% Checks
% % Get canonical rhs
[rhs_theta_c, rhs_hs_c, rhs_ctau_c] = canonical_rhs(BL, t, h, ctau, x);
% % Add Plasma terms to it
[rhs_theta_p, rhs_hs_p, rhs_ctau_p] = add_plasma_to_rhs(BL, t, h, ctau, x, rhs_theta_c, rhs_hs_c, rhs_ctau_c);
% % Transform rhs
[rhs_theta, rhs_h, rhs_ctau] = transform_rhs(BL, t, h, ctau, x, rhs_theta_p, rhs_hs_p, rhs_ctau_p);
% Done! Return!

% Now check transformation
[ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq);
[rhs_theta, rhs_h, rhs_ctau] = transformed_rhs(BL, t, h, ctau, x);


%% Integration 
% % Now explore integration
% Set initial point
y0(1) = t;
y0(2) = h;
y0(3) = ctau;

% Test function wrapper
[y_line] = BL.ode_fun_wrapper(x, y0);
% And handle function reference to wrapper
odefun = @BL.ode_fun_wrapper;
[y_line] = odefun(x, y0);

% Now test integration with simplest, low order, non-stiff ODE solver
% SOL = ode45(ODEFUN,[T0 TFINAL],Y0)
% Set integration range
x0    = xi_vector(1); xend  = xi_vector(end);

% Integrate
if fixed_step  == false
    % Variable Step Size
    sol = ode_solver(odefun,[x0 xend],y0);
else
    % Fixed Step Size
    sol = ode_solver(odefun,linspace(x0, xend, N_steps),y0);
end


%% Plotting 
% Extract some results
t_sol       = sol.y(1,:);
h_sol    = sol.y(2,:);
ctau_sol    = sol.y(3,:);
% And plot
figure(1)
plot(sol.x, h_sol)
% axis([0 100 0 4])
xlabel('x')
ylabel('H - Shape Factor')
figure(2)
plot(sol.x, re_theta( t_sol, u_inf_over_nu_inf ))
%axis([0 100 0 4])
xlabel('x')
ylabel('Re_theta')



%% Closure Function Summary
% % List of available functions
% % Definitions
% function [rt, rt_ue_over_nue ] = re_theta( t, ue_over_nue )
% function [ h, h_dstr, h_t] = H( dstr , t)
% function [ d , d_h, d_hh, d_t ] = delta_hh( h, hh, t)
% 
% % Closure Relations
% function [ cf, cf_hk, cf_rt, cf_msq  ] = cft( hk, rt, msq)
% function [ cf, cf_hk, cf_rt, cf_msq  ] = cft_rr( hk, rt, msq)
% function [ cq, cq_h , cq_hk , cq_hs , cq_us, cq_sv ] = cqt( h , hk , hs , us, sv )
% function [ cts, cts_h , cts_hk , cts_hs , cts_us, cts_sv ] = ctausuc( h , hk , hs , us, sv)
% function [ ctz, ctz_h , ctz_hk , ctz_hs , ctz_us ] = ctauzero( h , hk , hs , us)
% function [ di, di_hs, di_us, di_cf, di_st ] = dit( hs, us, cf, st)
% function [ di, di_hs, di_us, di_cf, di_st , di_sv] = dits( hs, us, cf, st, sv)
% function [ hc, hc_hk, hc_msq ] = hct( hk, msq )
% function [ hh, hh_hk ] = hh1cal( hk )
% function [ hk, hk_h, hk_msq  ] = hkin( h, msq)
% function [ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq)
% function [ us, us_h, us_hk, us_hs] = usg(h, hk, hs)
