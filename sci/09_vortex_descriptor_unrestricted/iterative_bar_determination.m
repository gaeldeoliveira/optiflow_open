% Study Script for determining bar variables from hat variables and tilde field


% % Load reference data for sandboxing and check consistency
% Load a converged solution for which we know all variables and fields: hat, tilde and
% bar
load final_state_alfa18

% Final Displayed State
% n = 10000    x = 1.5    CFL = 1.4041    PE = 5.3028
%      h_bar = 1.4             h_hat = 1.3341
%     hs_bar = 1.7444         hs_hat = 1.7766
%   dstr_bar = 0.010227      dstr_hat = 0.012475      dstr_tilde = 0.0022485
%  theta_bar = 0.007305     theta_hat = 0.0093511    theta_tilde = 0.0020461
% delta3_bar = 0.012743    delta3_hat = 0.016613    delta3_tilde = 0.0038706

% Test computation of shape factors
[hk_bar, hk_tilde, hk_hat] = hk_direct_external(   SMI , state.hk, state.rt, state.ue, nu_inf, state.u_tilde); %#ok<ASGLU>
%     hk_bar = 1.4000          hk_hat = 1.3341          hk_tilde = 1.0989

% Test joint computation of shape factors and rt numbers
[ hk_hat, rt_hat]         = hk_bar_rt_bar_external(SMI , state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
%    hk_bar = 1.4              hk_hat = 1.3341          hk_tilde = 1.0989
%    rt_bar = 500              rt_hat = 640.0452        rt_tilde = 140.0452
% Given
% state.rt = 500    and   state.hk = 1.4000



% % Proceed to make iterative method
% 
%       Purpose is to determine hk_bar and rt_bar (as if we didn't know them!)
%               from h_hat, rt_hat, nu_inf, ue and u_tilde
%       This problem is about finding the inverse of a function. We know
%       the direct function, because we have (or can construct easily):
%           a function for evaluating hk_hat and rt_hat 
%                from h_bar, rt_bar, nu_inf, ue and u_tilde
%       We will attempt the inversion with an iterative procedure. We will
%       choose a gradient 

% Settle reference hat values (inputs)
% For hk_hat
hk_hat_ref = hk_hat;
% For rt_hat
rt_hat_ref = rt_hat;


% Start with initial guesses
% For hk_bar
hk_bar0 = 1.45%; hk_hat;
% For rt_bar
rt_bar0 = 500 %; rt_hat;

% Set initial guess as current guess
hk_bar_i = hk_bar0;
rt_bar_i = rt_bar0;

% Beginning of future loop 

% Compute hk_hat and rt_hat values for current hk_bar and rt_bar guess
[ hk_hat_i, rt_hat_i]         = hk_bar_rt_bar_external(SMI , hk_bar_i , rt_bar_i , state.ue, nu_inf, state.u_tilde);

% Compute difference (to ref.)
delta_hk_hat = hk_hat_i - hk_hat_ref;
delta_rt_hat = rt_hat_i - rt_hat_ref;

% Compute violation

% Display violation
disp(['delta_hk_hat = ' num2str(   delta_hk_hat)]);
disp(['delta_rt_hat = ' num2str(   delta_rt_hat)]);

% Make new guess
hk_bar_i = hk_bar_i - delta_hk_hat ;
rt_bar_i = rt_bar_i - delta_rt_hat ;

disp(['hk_bar_i = ' num2str(hk_bar_i)]);
disp(['rt_bar_i = ' num2str(rt_bar_i)]);







