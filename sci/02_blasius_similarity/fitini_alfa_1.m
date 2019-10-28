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
N_steps     = 2000;                             % Set number of steps (only applicable for fixed_step mode
% fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps (always true for now, choice will be implemented later!)


%% Blasius Equation Definition
% Define Plasma Strenght
lambda = 0.1;
t_p    = 2;

% Define Boundary Condition Points
eta_A       = 0;
eta_B       = 20;
% Define Integration rage
eta_vector  = linspace(eta_A, eta_B, N_steps);

% Generate initial conditions structure
solinit = bvpinit(eta_vector,@odeinit_blasius);

% Solve (and pass lambda to fitini, la gentille équation!
sol = bvp_solver(@(eta, y) odefun_fitini(eta, y, lambda, t_p), @resBC_blasius , solinit);

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
 
 
legend('g - blasius function'       , ...
       'dg - velocity profile'      , ...
       'd2g - shear stress profile' , ...
       'Initial Guesses');

grid on
axis([0 5 -0.2 1.2])    % Manual Axis setting seems fairly universal!

%% Write

% Generate Speed Handle Function
u = @(eta) deval(sol, eta, 2);

% Generate Boundary Layer Parameters
dstr        = integral( @(eta)          1 - u(eta)     , eta_A, eta_B);
theta       = integral( @(eta) u(eta).*(1 - u(eta))    , eta_A, eta_B);
theta_str   = integral( @(eta) u(eta).*(1 - u(eta).^2) , eta_A, eta_B);
h           = dstr ./ theta;
h_str       = theta_str ./ theta;
cf_coef_x    = sol.y(3,1) * 2;
cf_coef_t    = sol.y(3,1) * 2 * theta;

% Display
disp(' 2*f_prime |  dstr    |   theta   |    h12   |    h_str    ');
disp(num2str([cf_coef_x , dstr , theta , h  , h_str]));




% % Auxiliary Functions (the mat4bvp example script uses nested functions,
% but I am afraid to make variable scope confusions with nested functions,
% so we will stay like this for now!)  

% % Some Results

% lambda = 0
%  2*f_prime |  dstr    |   theta   |    h12   |    h_str    
% 0.66412      1.7208      0.6641      2.5911      1.5726
% not hardcoded anymore!
% 0.66412      1.7208      0.6641      2.5911      1.5726

% lambda = -0.6
% b = 0.5 (t_p, but hardcoded!)
%  2*f_prime |  dstr    |   theta   |    h12   |    h_str    
% -0.33559      3.0869      1.0609      2.9097      1.5738
% not hardcoded anymore!
% -0.33559      3.0869      1.0609      2.9097      1.5738
% t_p = 1.0  (not hardcoded!)
% -0.3377      3.4338      1.0555      3.2533      1.5777

% lambda = 0.6
% b = 0.5 (t_p, but hardcoded!)
%  2*f_prime |  dstr    |   theta   |    h12   |    h_str    
% 4.4369     0.47507      0.2299      2.0664      1.6313
% not hardcoded anymore!
% 4.4369     0.47507      0.2299      2.0664      1.6313

