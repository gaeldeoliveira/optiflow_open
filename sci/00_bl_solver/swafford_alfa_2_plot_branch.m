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

% Clear Stuff
clear all; clc;


% Add some useful paths
addpath 0_closure_relations
addpath 1_definitions



% % % Profile Generation
% Step Numbers Follow Table 1 from supracited paper




% % Step 1
% Set Inputs

hk   = 3.5;                 % Shape Factor
rt  = 5000;                 % Reynolds Theta (Ue * theta / nue)
msq = 0;                    % Mach Number (for data consistency only!)

y_over_theta_vector = linspace(0,10,100);
u_over_ue_vector = zeros(size(y_over_theta_vector));

for n_y_over_theta = 1:length(y_over_theta_vector)
% % Step 2
% Calculate Skin Friction
[ cf, cf_hk, cf_rt, cf_msq  ] = cft( hk, rt, msq);
% Calculate ue_p = ue+
ue_p = sqrt(2 ./ abs(cf));

% % Step 3
% Compute Sign of Skin Friction
% s = cf / abs(cf)
s = sign(cf);

% % Step 4
% Compute u over ue (2)
u_ue_2 = (1 ./ 1.95) * ( atanh((8.5-hk)/7.5) - 0.364 );

% % Step 5
% Compute u over ue (5)
u_ue_5 = 0.155 + 0.795 * sech(0.51 * (hk-1.95));

% % Step 6
% Compute g(2)
g_2 = (u_ue_2 - (s ./ (0.09 * ue_p)) .* atan(0.18 * rt./ue_p)) ...
            ./ (1 - (s * pi()) ./ (0.18*ue_p) );
                
% % Step 7
% Compute g(5)
g_5 = (u_ue_5 - (s ./ (0.09 * ue_p)) .* atan(0.45 * rt./ue_p)) ...
            ./ (1 - (s * pi()) ./ (0.18*ue_p) );
        
% % Step 8
% Compute b (recall that in matlab log == ln)
b = log(atanh(g_2).^2 ./ atanh(g_5).^2) ./ log(2/5);

% % Step 9
% Compute a
a = atanh(g_2.^2) / (2.^b);

% % Step 10
% Evaluate profile!
y_over_theta = y_over_theta_vector(n_y_over_theta);
y_p = (rt / ue_p) .* y_over_theta;

u_p = (s ./ 0.09) * atan(0.09*y_p) + (ue_p - s*pi()/0.18) .* sqrt(tanh(a .* y_over_theta.^b));

u_over_ue = u_p / ue_p;
u_over_ue_vector(n_y_over_theta) = u_over_ue;

end

%% Plot

plot(u_over_ue_vector, y_over_theta_vector)
grid on
xlabel('U over Ue')
ylabel('y over theta')

%% Make derivative fit
% Make fit
pp = pchip(u_over_ue_vector, y_over_theta_vector);
% Take derivative
pp_du = ppdiff(pp);


% Set fitting properties
n_sources = 5;
y_top = 10;

% Fitting centers
y_origin_vector = ((1:n_sources)-0.5) / n_sources * y_top;
% Vorticity Integration Bounds
lower_integration_bounds = ((1:n_sources)-1) / n_sources * y_top;
upper_integration_bounds = ((1:n_sources)  ) / n_sources * y_top;




%% Legacy
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




