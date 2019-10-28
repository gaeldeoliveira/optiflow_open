% Study Script for determining bar variables from hat variables and tilde field


% % Load reference data for sandboxing and check consistency
% Load a converged solution for which we know all variables and fields: hat, tilde and
% bar

%load final_state_alfa18

% Final Displayed State
% n = 10000    x = 1.5    CFL = 1.4041    PE = 5.3028
%      h_bar = 1.4             h_hat = 1.3341
%     hs_bar = 1.7444         hs_hat = 1.7766
%   dstr_bar = 0.010227      dstr_hat = 0.012475      dstr_tilde = 0.0022485
%  theta_bar = 0.007305     theta_hat = 0.0093511    theta_tilde = 0.0020461
% delta3_bar = 0.012743    delta3_hat = 0.016613    delta3_tilde = 0.0038706

% Test computation of dstr and theta
[ dstr_bar,  dstr_tilde,  dstr_hat] = SMI.dstr_direct( state.hk, state.rt, state.ue, nu_inf, state.u_tilde);
[theta_bar, theta_tilde, theta_hat] = SMI.theta_direct(state.hk, state.rt, state.ue, nu_inf, state.u_tilde);

% Display Values
disp(['  dstr_bar = ' num2str( dstr_bar)  '      dstr_hat = ' num2str(  dstr_hat) '      dstr_tilde = ' num2str(  dstr_tilde)]);
disp([' theta_bar = ' num2str( theta_bar) '     theta_hat = ' num2str( theta_hat) '     theta_tilde = ' num2str( theta_tilde)]);

%  dstr_bar = 0.010227      dstr_hat = 0.012475      dstr_tilde = 0.0022485
% theta_bar = 0.007305     theta_hat = 0.0093511     theta_tilde = 0.0020461
% Given
% state.rt = 500    and   state.hk = 1.4000



% % Proceed to make iterative method
% 
%       Purpose is to determine dstr_bar and theta_bar (as if we didn't know them!)
%               from h_hat, rt_hat, nu_inf, ue and u_tilde
%       This problem is about finding the inverse of a function. We know
%       the direct function, because we have (or can construct easily):
%           a function for evaluating dstr_hat and theta_hat 
%                from h_bar, rt_bar, nu_inf, ue and u_tilde
%       We will attempt the inversion with an iterative procedure. We will
%       choose a gradient 

% Settle reference hat values (inputs)
% For dstr_hat
dstr_hat_ref  = dstr_hat;
% For theta_hat
theta_hat_ref = theta_hat;


% Start with initial guesses
% For dstr_bar
dstr_bar0  = dstr_hat;
% For theta_bar
theta_bar0 = theta_hat;

% Set initial guess as current guess
dstr_bar_i  = dstr_bar0;
theta_bar_i = theta_bar0;

% Beginning of future loop 
for n=1:10
disp(['New Iteration   n = ' num2str( n)]);
% Compute hk_bar_i and rt_bar_i values for current dstr_bar and theta_bar guess
hk_bar_i = dstr_bar_i / theta_bar_i;
rt_bar_i = theta_bar_i * state.ue / nu_inf;

% Display
disp(['hk_bar_i = ' num2str( hk_bar_i)]);
disp(['rt_bar_i = ' num2str( rt_bar_i)]);

% Compute dstr_hat and dstr_hat values for current hk_bar and rt_bar guess
[ dstr_bar_i,  dstr_tilde_i,  dstr_hat_i] = SMI.dstr_direct( hk_bar_i , rt_bar_i , state.ue, nu_inf, state.u_tilde);
[theta_bar_i, theta_tilde_i, theta_hat_i] = SMI.theta_direct(hk_bar_i , rt_bar_i , state.ue, nu_inf, state.u_tilde);

% Display Values
disp(['  dstr_bar = ' num2str( dstr_bar)  '      dstr_hat = ' num2str(  dstr_hat) '      dstr_tilde = ' num2str(  dstr_tilde)]);
disp([' theta_bar = ' num2str( theta_bar) '     theta_hat = ' num2str( theta_hat) '     theta_tilde = ' num2str( theta_tilde)]);

% Compute difference (to ref.)
delta_dstr_hat  =  dstr_hat_i -  dstr_hat_ref;
delta_theta_hat = theta_hat_i - theta_hat_ref;

% Compute violation

% Display violation
disp(['delta_dstr_hat = ' num2str(   delta_dstr_hat)]);
disp(['delta_theta_hat = ' num2str(   delta_theta_hat)]);

% Make new guess
dstr_bar_i  =  dstr_bar_i -  delta_dstr_hat ;
theta_bar_i = theta_bar_i - delta_theta_hat ;

% Display New Guess
disp([' dstr_bar_i = ' num2str( dstr_bar_i)]);
disp(['theta_bar_i = ' num2str(theta_bar_i)]);





end

