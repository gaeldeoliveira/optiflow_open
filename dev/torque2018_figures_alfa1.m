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
CS0 = open('fig/run-over-database-alfa-15/mat/run_over_database_alfa15_20171123T100928.mat');
% Case with post-processing (seems to be the right one!)
CS1 = open('fig/run-over-database-alfa-15/mat/run_over_database_alfa15_20171123T105445.mat');
%
CS1.xgrad

% Now prepare plot of Cf closure relation
msq      = 0;
re_theta = 1000;
hk_range = linspace(1,10);

% Compute Prelearn Reference Skin Friction Relation 
[ cf_reference] = cft(    hk_range, re_theta , msq);
% Compute Prelearn Increased Skin Friction Relation
[ cf_increased] = cft_rr( hk_range, re_theta , msq);
% Compute Postlearn Increased Skin Friction Relation
[ cf_postlearn] = cft_ml( hk_range, re_theta , msq, CS1.SPR);

% Plot Comparison of prelearn Cf relations
plot(hk_range, cf_reference, hk_range, cf_increased, hk_range, cf_postlearn); grid on;

% Plot 
