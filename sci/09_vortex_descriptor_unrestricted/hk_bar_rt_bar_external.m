function [ hk_hat, rt_hat] = hk_bar_rt_bar_external( SMI, hk_bar, rt_bar, ue, nue, u_tilde)
%H_DIRECT_EXTERNAL computes the h shape factor and the rt number for the 
% hat field, given a ShearMixedIntegrator (SMI) object, the shear field 
% state (hk_bar, rt_bar), edge flow conditions (ue, nue) and the u_tilde 
% field state (a 2-d mesh consistent with SMI.CM CrossflowMesh object)
%
% The function resorts to the standard routines of the SMI object :
%       SMI.dstr_direct(...)
%       SMI.theta_direct(...)
% This approach is far from efficient, but it does the job with minimal
% effort, which is what are looking for at this stage!
%

% Compute the displacement thicknesses (dstr) for the three fields (bar, tilde, hat)
[dstr_bar, dstr_tilde, dstr_hat] = dstr_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde);

% Compute the momentum thicknesses (theta) for the three fields (bar, tilde, hat)
[theta_bar, theta_tilde, theta_hat] = theta_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde);

% Compute the shape factors
hk_bar     =   dstr_bar ./ theta_bar  ;     % bar   field (should match inputs exactly!)
hk_tilde   = dstr_tilde ./ theta_tilde;     % tilde field (this is the one that only depends on integrations!)
hk_hat     =   dstr_hat ./ theta_hat  ;     % hat   field (depends on reconstructed inputs (bar component) and integrations (tilde component))
% Care should be taken to settle the differences in meaning of tilde
% integral variables:
%    -> For the thickensses (dstr, theta, d3) the tilde component is
%       the difference between the hat and bar components
%    -> For the shape factors, the tilde component is a ratio between tilde
%       thicknesses. The meaning of which is not obvious (not its use, or 
%       even whether it has any, in fact!)

% Compute the rt numbers 
rt_bar     = theta_bar   * ue / nue;
rt_tilde   = theta_tilde * ue / nue;
rt_hat     = theta_hat   * ue / nue;

% Display current computation

disp(['    hk_bar = ' num2str(   hk_bar)  '        hk_hat = ' num2str(    hk_hat) '        hk_tilde = ' num2str(    hk_tilde)]);
disp(['    rt_bar = ' num2str(   rt_bar)  '        rt_hat = ' num2str(    rt_hat) '        rt_tilde = ' num2str(    rt_tilde)]);


end

