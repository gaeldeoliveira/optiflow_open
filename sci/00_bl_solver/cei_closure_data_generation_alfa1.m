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
%       Energy Interaction Coefficient Closure Data Generator
%           following:
%               Modelling the effect of DBD plasma actuators on boundary
%               layer development
%               Internal Report, Gael de Oliveira, Ricardo Pereira
%               September 8, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
% Clear Stuff
close all; clear all; clc;

% Add some useful paths
addpath 0_closure_relations
addpath 1_definitions
addpath 2_helper_functions

%% Inputs

% % Data Generation Sets
% Set of shape factors (hk) 
hk_range = linspace(1.4, 8, 3); 
% Set of Momentum Thickness Reynolds (Re_theta = rt)
rt_range = linspace(500, 5000, 3); 
% Set of Momentum Scaled Plasma Force Field Coefficients (t_theta_p)
t_theta_p_range = linspace(1.4, 8, 3); 

% % Numerics parameters
N_integration_points = 1000;            % Integration seems very stable with naive middle point riemann scheme, for 100 points or more (1000 is still a good measure, specially for low shape factors)


%% Object Instanciation
% Construct Swafford Profile Object
SP = swafford_profile();

% Construct Plasma Descriptor Object (used for wy_function only)
PD = plasma_descriptor();

% Construct cei integrand function
cei_integrand = @(yt, b) SP.evaluate_profile(yt) .* PD.wy_function(yt, b) ./ b;

%% Data Allocation
% Allocate case arrays for each study dimension
N_hk        = length(hk_range);
N_rt        = length(rt_range);
N_t_theta_p = length(t_theta_p);
