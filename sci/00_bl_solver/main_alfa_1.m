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


% List of available functions
% Definitions
function [rt, rt_ue_over_nue ] = re_theta( t, ue_over_nue )
function [ h, h_dstr, h_t] = H( dstr , t)

% Closure Relations
function [ cf, cf_hk, cf_rt, cf_msq  ] = cft( hk, rt, msq)
function [ cf, cf_hk, cf_rt, cf_msq  ] = cft_rr( hk, rt, msq)
function [ cq, cq_h , cq_hk , cq_hs , cq_us, cq_sv ] = cqt( h , hk , hs , us, sv )
function [ cts, cts_h , cts_hk , cts_hs , cts_us, cts_sv ] = ctausuc( h , hk , hs , us, sv)
function [ ctz, ctz_h , ctz_hk , ctz_hs , ctz_us ] = ctauzero( h , hk , hs , us)
function [ di, di_hs, di_us, di_cf, di_st ] = dit( hs, us, cf, st)
function [ di, di_hs, di_us, di_cf, di_st , di_sv] = dits( hs, us, cf, st, sv)
function [ hc, hc_hk, hc_msq ] = hct( hk, msq )
function [ hh, hh_hk ] = hh1cal( hk )
function [ hk, hk_h, hk_msq  ] = hkin( h, msq)
function [ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq)
function [ us, us_h, us_hk, us_hs] = usg(h, hk, hs)


% Draft RHS equation evalutation
% Set inputs
re_c        = 3e6
ue          = 1.0;
due_dxi     = 0;                % Convention total derivatives get a d in there!
msq         = 0;
ue_over_nue = re_c * 0.8;       % Take 80% chord stance


% Receive particular request
t       = 0.002;
dstr    = 0.003;
ctau    = 0.06*0.06;

% % % Start Calculation
% Shape Factor
[ h, h_dstr, h_t] = H( dstr , t);                           % Looks ok
% Kinematic Shape Factor
[ hk, hk_h, hk_msq  ] = hkin( h, msq);                      % Looks ok

% Reynolds Theta
[rt, rt_ue_over_nue ] = re_theta( t, ue_over_nue );         % Looks ok

% Now calculate Cf
[ cf, cf_hk, cf_rt, cf_msq  ] = cft_rr( hk, rt, msq);       % Looks ok

% Momentum Equation RHS Evaluation Draft
rhs_momentum = cf/2 + (h + 2 - msq.^2) .* (t ./ ue) .*  due_dxi;

% Now energy shape factor (H_star)
[ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq)


