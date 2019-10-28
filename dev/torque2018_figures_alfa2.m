%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   LearneRfoil Development Framework
%
%   Author:
%               Gael de Oliveira November/December 2017   
%
%   File Description:
%               torque2018 processing script for results from
%               run_over_database_alfa15 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add Necessary Paths
addpath src
addpath user_src
addpath dev
addpath dev/0_closure_relations/

% Case without post-processing (not fully converged, not sure)
%CS0 = open('fig/run-over-database-alfa-15/mat/run_over_database_alfa15_20171123T100928.mat');
% Case with post-processing (seems to be the right one!)
%CS1 = open('fig/run-over-database-alfa-15/mat/run_over_database_alfa15_20171123T105445.mat');
CS1 = open('fig/run-over-database-alfa-17/mat/run_over_database_alfa17_intermediate_processing_20171204T045816.mat');

% Insert Experimental Data
rt_vinuesa = 2255  ; %#ok<NASGU>
hk_vinuesa = 2.03  ; %#ok<NASGU>
cf_vinuesa = 0.0012; %#ok<NASGU>

rt_vinuesa = 1722  ;
hk_vinuesa = 1.74  ;
cf_vinuesa = 0.0024;

% Now prepare plot of Cf closure relation
msq      = 0;
re_theta = rt_vinuesa;
hk_range = linspace(1,5);

% Compute Prelearn Reference Skin Friction Relation 
[ cf_reference] = cft(    hk_range, re_theta , msq);
% Compute Prelearn Increased Skin Friction Relation
[ cf_increased] = cft_rr( hk_range, re_theta , msq);
% Compute Postlearn Increased Skin Friction Relation
[ cf_postlearn] = cft_ml( hk_range, re_theta , msq, CS1.SPR);

% Plot Comparison of prelearn Cf relations
figure(87)
subplot(221)
%plot(hk_range, cf_reference, hk_range, cf_increased, hk_range, cf_postlearn); grid on;
plot(hk_range, cf_reference * 1e4, hk_range, cf_postlearn* 1e4, hk_vinuesa, cf_vinuesa* 1e4, 'kx'); grid on;
legend('Previous Closure', 'Learned Closure', 'DNS [Vinuesa et al.]');
xlabel('H_k - Kinematic Shape Factor');
ylabel('C_f x 10^4');
title(['Skin Friction (Re_{\theta} = ' num2str(re_theta) , ')']);
axis([1 5 -20 80])

% % Now go on to explore H* correlation

% Compute Prelearn Reference H* Skin Friction Relation 
[ hs_reference] = hst(    hk_range, re_theta , msq);
% Compute Prelearn Increased H* Skin Friction Relation
[ hs_increased] = hstivw( hk_range, re_theta , msq);
% Compute Postlearn Increased H* Skin Friction Relation
[ hs_postlearn] = hst_ml( hk_range, re_theta , msq, CS1.SPR);

% Plot Comparison of prelearn Cf relations
% figure(85)
subplot(222)
%plot(hk_range, hs_reference, hk_range, hs_increased, hk_range, hs_postlearn); grid on;
plot(hk_range, hs_reference, hk_range, hs_postlearn); grid on;
legend('Previous Closure', 'Learned Closure');
xlabel('H_k - Kinematic Shape Factor');
ylabel('H^{*} - Skin Friction Coefficient');
title(['Energy Shape Factor (Re_{\theta} = ' , num2str(re_theta), ')']);
axis([1 5 1.5 2])


print(['./closure_relations_improved' , datestr(now(), 30)], '-dpdf')




% Plot 
