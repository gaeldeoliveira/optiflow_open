%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Integral Boundary Layer Integrator (ODE)
%           Plasma Development Tool
%
%       August 2014
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

% Plot H_star vs H as some Re_theta
msq = 0 ; rt = 300 ; hk = linspace(1, 8, 400);
hs = zeros(size( hk));

for n_hk = 1:length(hk)
    hs(n_hk)= hst( hk(n_hk), rt, msq);
end

plot(hk, hs)


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


% Draft RHS equation evalutation
% Set inputs
xi_vector   = linspace(0, 10, 100);
ue_vector   = ones(size(xi_vector));
msq_vector  = zeros(size(xi_vector));

% Define nue
re_c        = 3e6;
ue_over_nue = re_c;       % Take 80% chord stance
nue         = 1 ./ ue_over_nue;

% Create Boundary Layer Object
BL = boundary_layer();
% Set nue and forcing terms (order matters!)
BL.nue = nue;
BL.set_forcing_terms(xi_vector, ue_vector, msq_vector)

% Define particular RHS request for testing
t       = 0.002;
dstr    = 0.003;
ctau    = 0.0015;         % Max shear stress coefficient
x       = 0.0;

[rhs_theta, rhs_hs, rhs_ctau] = canonical_rhs(BL, t, dstr, ctau, x);
[rhs_theta, rhs_dstr, rhs_ctau] = transformed_rhs(BL, t, dstr, ctau, x)

% % Now explore integration
% Set initial point
y0(1) = t ;
y0(2) = dstr;
y0(3) = ctau;

% Test function wrapper
[y_line] = BL.ode_fun_wrapper(x, y0)
% And handle function reference to wrapper
odefun = @BL.ode_fun_wrapper;
[y_line] = odefun(x, y0)

% Now test integration with simplest, low order, non-stiff ODE solver
% SOL = ode45(ODEFUN,[T0 TFINAL],Y0)
% Set integration range
x0    = 0; xend  = 1;
% Integrate
sol = ode45(odefun,[x0 xend],y0);

% Extract some results
t_sol       = sol.y(1,:);
dstr_sol    = sol.y(2,:);
ctau_sol    = sol.y(3,:);

% And Post process them (should tend towards h = 1.4!)
h_sol       = dstr_sol ./ t_sol;




plot(sol.x, h_sol)