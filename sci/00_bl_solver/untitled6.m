% Script to show us - cf inconsistency near separation


%% Inputs
ue_over_nue = 1e6;                  % Reynolds Number Proxy
t           = 0.02;                 % Theta - Displacement Thickness
msq         = 0;                    % Mach number (0 for incompressible)
h = linspace(1.2, 5, 100);          % Shape factor range

%% Compute curves
% Reynolds Theta
[rt, rt_ue_over_nue ] = re_theta( t, ue_over_nue );




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
