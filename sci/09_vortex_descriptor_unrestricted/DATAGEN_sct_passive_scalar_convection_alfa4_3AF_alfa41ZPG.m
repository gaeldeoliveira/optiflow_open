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
%TEMPORARY
%Regular      : state.y_v     = h_VG;             % [m /s ] distance of filament core center to wall
%Current trial: state.y_v     = 0.6 * h_VG;             % [m /s ] distance of filament core center to wall


%% Preparation
clear all; close all;clc; 

% Add Paths
addpath closure_relations/
addpath definitions/
addpath intersections/


vg_case.rho_inf         = 1.225;                  % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.461e-5;               % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 15.16;                  % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                   % [adim.] Shape Factor
vg_case.rt0             = 2498;                   % [adim.] Reynolds Theta of Starting Flow

vg_case.D_VG            =  30.0e-3;                 % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  12.5e-3;                 % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  12.5e-3;                 % [m    ] chord of VG  (5cm)
vg_case.h_VG            =   5.0e-3;                 % [m    ] height of VG
vg_case.AoA_VG          = -18     ;                 % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    25    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.50 * vg_case.h_VG;   % Time between information events

vg_case.snapshot_period = 40;

%% Inputs



