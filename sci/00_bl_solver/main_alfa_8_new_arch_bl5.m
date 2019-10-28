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
close all; clear all; clc;


% Add some useful paths
addpath 0_closure_relations
addpath 1_definitions

%% Inputs
% Initial Conditions (ctau will be guessed later)
x       = 0.0;                          % Starting Stance
t       = 0.003;                        % Momentum Thickness
h       = 2.5;                          % Shape Factor
msq     = 0;                            % Mach Number

% Define Forcing Terms
xi_vector   = linspace(0, 100, 100);
ue_vector   = ones(size(xi_vector));
msq_vector  = zeros(size(xi_vector));

% Define viscosity (nue)
re_c        = 3e6;
ue_over_nue = re_c;                 % Take full chord (1) as reference
nue         = 1 ./ ue_over_nue;


%% Allocation and Instanciation
% Create Boundary Layer Object
BL = boundary_layer5();

% Set nue 
BL.nue = nue;

% Set forcing terms (order matters! must be done after nue)
BL.set_forcing_terms(xi_vector, ue_vector, msq_vector)

%% Build Consistent Initial Conditions
% And guess corresponding ctau from equilibrium ctau (even though flat plate is not in perfect equilibrium!)
hk      = hkin( h, msq);                    
rt      = re_theta( t, ue_over_nue );       % Looks ok
cf      = cft( hk, rt, msq);                % Looks ok
hs      = hst( hk, rt, msq);
us      = usg(h, hk, hs);
ctau    = ctauzero( h , hk , hs , us);      % Looks ok
cdi     = cdissipation( us, cf, ctau);      % Looks ok
% Cannonical RHS
[rhs_theta, rhs_hs, rhs_ctau] = canonical_rhs(BL, t, h, ctau, x);

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
[y_line] = BL.ode_fun_wrapper(x, y0)
% And handle function reference to wrapper
odefun = @BL.ode_fun_wrapper;
[y_line] = odefun(x, y0)

% Now test integration with simplest, low order, non-stiff ODE solver
% SOL = ode45(ODEFUN,[T0 TFINAL],Y0)
% Set integration range
x0    = xi_vector(1); xend  = xi_vector(end);
% Integrate
sol = ode45(odefun,[x0 xend],y0);


%% Plotting 
% Extract some results
t_sol       = sol.y(1,:);
h_sol    = sol.y(2,:);
ctau_sol    = sol.y(3,:);
% And plot
figure(1)
plot(sol.x, h_sol)
axis([0 100 0 4])
xlabel('x')
ylabel('H - Shape Factor')
figure(2)
plot(sol.x, re_theta( t_sol, ue_over_nue ))
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
