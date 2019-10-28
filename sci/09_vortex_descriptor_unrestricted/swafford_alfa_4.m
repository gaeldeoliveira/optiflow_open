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
%       Swafford velocity profile generator
%               following:
%                   Analytical Approximation of Two-Dimensional Separated
%                   Turbulent Boundary-Layer Velocity Profiles
%                   Swafford, T.W., AIAA Journal Vol.21, N.6, June 1983
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Script to generate plots:
%               Swafford Profile
%               Cei Integral Convergence Study
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clear Stuff
close all; clear all; clc;


% Add some useful paths
addpath 0_closure_relations


% Construct Swafford Profile
SP = swafford_profile();

% Set hk and rt


% Plot Profile
figure(1)
hk = 1.4; rt = 500;
SP.plot(hk, rt, 'b')
hold on
hk = 1.4; rt = 5000;
SP.plot(hk, rt, 'b--')
hk = 6.0; rt = 500;
SP.plot(hk, rt, 'r')
hk = 6.0; rt = 5000;
SP.plot(hk, rt, 'r--')
legend('H_k = 1.4     Re_{theta} = 500' , 'H_k = 1.4     Re_{theta} = 5000', ... 
       'H_k = 6       Re_{theta} = 500' , 'H_k = 6       Re_{theta} = 5000', ...
       'Location', 'NorthWest');
title('Swafford Profile');
xlabel('U/U_e') ; ylabel('y/delta_2')

set(gcf, 'PaperType', 'A5'); orient landscape
%print -dpdf -painters ./figures/swafford_profile.pdf

% Set thickness integration bounds
yt_min = 0;
yt_max = 100;

% Use Gauss Kronberg Quadrature for integrals
quad_scheme = @integral;

% Update back to ref conditions
hk = 1.4; rt = 500;
SP.plot(hk, rt, 'b')

% Define Delta Star Integrand
dstr_integrand  = @(yt) (1 - SP.evaluate_profile(yt));
theta_integrand = @(yt) SP.evaluate_profile(yt) .* (1 - SP.evaluate_profile(yt));

% Integrate 
% [dstr  edstr ] = quad_scheme(@(yt) dstr_integrand(yt) , yt_min, yt_max);
% [theta etheta]  = quad_scheme(@(yt) theta_integrand(yt) , yt_min, yt_max);


% Straight integration
y_integration = linspace(yt_min, yt_max, 10000000);
% Cosine distribution integration
% phi_integration = linspace(0, pi()/2, 10000000);
% y_integration   = [2*(1 - cos(phi_integration)) , 2 + phi_integration * 2 * yt_max / pi()];


dstr  = riemann_integral(@(yt) dstr_integrand(yt), y_integration);
theta  = riemann_integral(@(yt) theta_integrand(yt), y_integration);
[dstr , theta , dstr / theta , (dstr / theta - hk) / hk]


% Define integral 

wy_function = @(y, b) 0.5*pi() * sin( 0.5*pi()*(y/b + 1)) .* (y/b<1);;
b = 8;
cei_integrand = @(yt) SP.evaluate_profile(yt) .* wy_function(yt, b) ./ b;

N_points_log_10  = 1:7;
N_points         = 10.^(1:7);
cei_vector1      = zeros(size(N_points));
cei_vector2      = zeros(size(N_points));
cei_vector3      = zeros(size(N_points));
cei_vector4      = zeros(size(N_points));
% Test numerical stability 
for n=1:length(N_points)
    hk = 1.4; rt = 500; SP.update_hk_rt_pair(hk, rt);
    y_integration = linspace(0, b, N_points(n));
    cei_vector1(n) = riemann_integral(@(yt) cei_integrand(yt), y_integration);
    
    hk = 1.4; rt = 5000; SP.update_hk_rt_pair(hk, rt);
    y_integration = linspace(0, b, N_points(n));
    cei_vector2(n) = riemann_integral(@(yt) cei_integrand(yt), y_integration);
    
    hk = 6; rt = 500; SP.update_hk_rt_pair(hk, rt);
    y_integration = linspace(0, b, N_points(n));
    cei_vector3(n) = riemann_integral(@(yt) cei_integrand(yt), y_integration);
    
    hk = 6; rt = 5000; SP.update_hk_rt_pair(hk, rt);
    y_integration = linspace(0, b, N_points(n));
    cei_vector4(n) = riemann_integral(@(yt) cei_integrand(yt), y_integration);    
end

% Conclusions:
%   1000 equally spaced integration points is more than enough
%   The profile becomes imaginary for small (rt<1000 and hk<1.4) input
%   pairs. Routine to eliminate imaginary results is needed!

% Show plot for small Re
figure(2)
plot(N_points_log_10, cei_vector1 , N_points_log_10, cei_vector3, ...
     N_points_log_10, cei_vector3 , N_points_log_10, cei_vector4);
grid on
xlabel('log_{10}N_{points}');
ylabel('C_{ei} - Energy Interaction Coefficient')
legend('H_k = 1.4     Re_{theta} = 500' , 'H_k = 1.4     Re_{theta} = 5000', ... 
       'H_k = 6       Re_{theta} = 500' , 'H_k = 6       Re_{theta} = 5000');
axis([1 7 -.2 1])
title('Riemann Integration Convergence -- t_theta_p = 8')

set(gcf, 'PaperType', 'A5'); orient landscape
%print -dpdf ./figures/cei_integral_convergence.pdf
   

%plot(linspace(0,2,100), wy_function(linspace(0,2,100),1.5))




