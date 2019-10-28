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


%% Inputs

% Integration algorithm parameters
options.bvp_solver  = @bvp4c;                           % bvp4c and bvp5c are the two common choices
options.N_steps     = 2000;                             % Set number of steps (only applicable for fixed_step mode
% fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps (always true for now, choice will be implemented later!)

% Define Plasma Strenght
options.lambda = 0;
options.t_p    = 2;

% Define Boundary Condition Points (range of integration, eta_B ~= infty)
options.eta_A       = 0;
options.eta_B       = 20;

%% Solve Fitini
outputs = solver_fitini(options);

%% Display Solution
disp('| lambda  |   t_p   | 2*f_prime |  dstr    |   theta   |    h12   |    h_str    ');
disp([ '    ' ,  num2str([options.lambda options.t_p outputs.cf_coef_x , outputs.dstr , outputs.theta , outputs.h  , outputs.h_str])]);


%% Plot
plot(outputs.eta_vector , outputs.y_sol(1,:)  , ...
     outputs.eta_vector , outputs.y_sol(2,:)  , ...
     outputs.eta_vector , outputs.y_sol(3,:) );

hold on
 
plot(outputs.eta_vector , outputs.y_0(1,:)  , '--' , ...
     outputs.eta_vector , outputs.y_0(2,:)  , '--' ,  ...
     outputs.eta_vector , outputs.y_0(3,:)  , '--');
 
 
legend('g - blasius function'       , ...
       'dg - velocity profile'      , ...
       'd2g - shear stress profile' , ...
       'Initial Guesses');

grid on
axis([0 5 -0.2 1.2])    % Manual Axis setting seems fairly universal!

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


