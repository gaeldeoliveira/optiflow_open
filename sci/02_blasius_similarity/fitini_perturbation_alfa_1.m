%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Similarity Profile Solver (ODE)
%           Plasma Development Tool
%
%       October 2014, GNU-GPLv3 or later
%       Gael de Oliveira, Ricardo Pereira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all

% function blasius_alfa_1()

%% Inputs

% Integration algorithm parameters
bvp_solver  = @bvp4c;                           % bvp4c and bvp5c are the two common choices
N_steps     = 1000;                             % Set number of steps (only applicable for fixed_step mode
% fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps (always true for now, choice will be implemented later!)


%% Blasius Equation Definition
% Define Boundary Condition Points
eta_A       = 0;
eta_B       = 10;
% Define Integration rage
eta_vector  = linspace(eta_A, eta_B, N_steps);

% Generate initial conditions structure
solinit = bvpinit(eta_vector,@odeinit_blasius);

% Solve
sol = bvp_solver(@odefun_fitini , @resBC_blasius_fitini_perturbation , solinit);

% end

% Now retrieve solution
y_sol = deval(sol,eta_vector);
y_0   = solinit.y;

%% Plot
plot(eta_vector , y_sol(1,:)  , ...
     eta_vector , y_sol(2,:)  , ...
     eta_vector , y_sol(3,:) );

hold on
 
plot(eta_vector , y_0(1,:)  , '--' , ...
     eta_vector , y_0(2,:)  , '--' ,  ...
     eta_vector , y_0(3,:)  , '--');
 
 
legend('g - blasius function'  , ...
       'dg - velocity profile' , ...
       'd2g - shear stress profile');

grid on


% % Auxiliary Functions (the mat4bvp example script uses nested functions,
% but I am afraid to make variable scope confusions with nested functions,
% so we will stay like this for now!)  

