%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Kirikou -   A simple 2d Vorticity Equation Solver for                 %
&               actuator disk flows                                       %
%                                                                         %
%   Author  :   Gael de Oliveira                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clean Environment (a bit brutally, for sure!)
clear all; close all; clc; %#ok<CLALL>


% % Inputs are adimensionalized to references:
%       u_inf = 1     m/s       Freestream Speed
%       rho   = 1.225 kg/m3     Fluid Density (incomp.)
%       d     = 1     m         Diameter
Ct                  =  8/9;     % [-- ] Actuator Force Coefficient (C_F_a)
x_v                 = -0.2;     % [-- ] x-stance of lifting vortex pair (symmetry on x axis)
r_v                 =  0.4;     % [-- ] distance of lifting vortex pair to x axis (symmetry on x axis)
gamma_v             =  0.0;     % [-- ] strenght of each vortex in pair

% % Create Solver Object
DSD = kirikou_single_actuator_solver(Ct, x_v, r_v, gamma_v);

DSD.information_period = 100; 
DSD.n_stances          = 20;
DSD.relax              = 0.001;                 % Stretching Relaxation
DSD.dt                 = 0.001;                 % Shape (Virtual) Time Step

% % Solve!
DSD.preprocess_run_postprocess();

% % Enjoy results!
DSD.plot_velocity_field_with_streamlines();

% % Explore the fields of the DSD object, with a particular focus on the
% DSD.RES structure, which summarizes case results!
% Typical DSD.RES values (for example inputs):
%             r_a: 0.5000               % Actuator Radius
%           phi_a: -0.4444              % Actuator Force (Loading) Density
%             x_v: -0.2000              % x-stance of lifting vortices
%             r_v: 0.4000               % y-stance of lifting vortices
%         gamma_v: 0.4000               % circulation strenght of lifting vortices
%           C_F_a: 0.8889               % Actuator Force Coefficient
%           C_F_b: 0.2088               % Vortices Streamwise Force Coefficient (computed with generalized Kutta-Joukowski-Lagally theorem using numerical velocity)
%          Cp_num: -0.7348              % Numerical Cp, computed by integrating velocity field over actuator
%         Cp_theo: -0.7318              % Theoretical Cp, computed from C_F_b using the de Vries power coefficient law (as in my Torque 2016)
%            e_Cp: 0.0041               % Absolute Difference between Numerical Cp values
%          ua_num: 0.8266               % Average velocity on actuator plane (numerical)
%         ua_theo: 0.8233               % Average velocity on actuator plane (theorethical, given C_F_a and C_F_b)
%          a1_num: 0.1734               % Induction Factor on actuator plane (numerical)
%         a1_theo: 0.1767               % Induction Factor on actuator plane (theoretical)
%     RES_stretch: 9.9787e-13           % Solution Residual Filament Vorticity Evolution (RMS)
%       RES_shape: 8.8902e-13           % Solution Residual Wake Flow Tangency Condition (RMS)
%              VS: [1x1 constant_strenght_vortex_segment_2d]
