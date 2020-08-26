%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Dogoro -    A simple specialized 2d Vorticity Equation Solver for %
%                   Actuator Disk Flows                                   %
%                                                                         %
%       Date    :   June 2014 to October 2015                             %
%       Authors :   Gael de Oliveira                                      %
%                   Ricardo Pereira                                       %
%       License :   GNU GPLv2 or later                                    %
%                                                                         %
%       If you use this work for your research, you can cite this article:%
%                                                                         %
%       Oliveira, G.L, Pereira, R.B. & Van Bussel, G.B., "Can we learn    %
%       something by placing a wind turbine behind a fan?", Nature Energy,%
%       xx/2016                                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       File    :   main_cpt_illustration_alfa3.m                         %
%       Title   :   Make Plots of the Analytical Momentum Model           %
%                                                                         %
%       Purpose :   This script is used to plot the asymptotic momentum   %
%                   model predictions of the power coefficient and optimal%
%                   operation regimes.                                    %
%                   This code is self-contained in a single file thanks   %
%                   to the use of inline functions. Many functions are    %
%                   declared twice (in theo and num version) to allow for %
%                   easy verification of the algebraic manipulations and  %
%                   avoid/correct syntax errors!                          %
%       Note    :   Strictly speaking, this is not part of the solver.    %
%                   In fact, it illustrates the model that we want to     %
%                   validate with the vorticity solver.                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean Environment (a bit brutally, for sure!)
clear all; close all; clc; %#ok<CLALL>

% Inputs
% Environment
u_inf = 1;       % Unperturbed Wind Speed
rho   = 1;       % Density
% Machine Parameters
a2_in = 1/3;                                        % 0 < a2_in < 0.5
b_in  = 1.5;                                        % b>0 (=1 means no actuator at front, >1 means energy addition on front!)
Dx    = 2.0;
% Rotor Forces (prescribed, generated from a1_in and b_in parameters!)
f1  = - 0.5 * rho * u_inf.^2 * (b_in^2 - 1);
f2  = - 0.5 * rho * u_inf.^2 * (b_in^2 * ((1-2*a2_in)^2-1));

% Instanciate Solver Object
DSD = dogoro_double_actuator_solver(f1, f2, Dx);

% Make sure environment variables are consistent
DSD.u_inf = u_inf;
DSD.rho   = rho;

% Set Convergence Parameters
DSD.information_period = 400;
DSD.N_iterations       = 40000;
DSD.relax_shape        = 0.005;
DSD.relax_stretch      = 0.005;

% Solve
DSD.preprocess_run_postprocess();

% Plot and enjoy!
DSD.plot_velocity_field_with_streamlines();
