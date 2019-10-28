% Clean environment
clear all; close all; clc;

% System Work necessary to start
fs = filesep();                 % Folder separator is OS dependent
addpath([cd fs 'src']);         % Add optimizer source code folder
addpath([cd fs 'user_src']);    % Add optimizer user sources (cost functions, etc) folder
addpath([cd fs 'gui']);         % Add gui source code folder
addpath([cd fs 'dev']);         % Add devsource code folder
warning('off','MATLAB:DELETE:FileNotFound');
warning('off','MATLAB:MKDIR:DirectoryExists');


% Get a database
loaded_datase = load('./data/experimental/EDB18dfdd03-e00b-4717-8af4-6ba3892db2e5.mat');
% Make it present
EDB = loaded_datase.EDB;
% Restrict to a subset (for faster prototyping)
EC_cell_full = EDB.EC_cell;
%EDB.EC_cell = EC_cell_full(1:10);
EDB.EC_cell = EC_cell_full(13:60); % Very thin airfoils (first ones, 0006, 0009, 1408 and 1410 lead to issues (convergence?), probably due to LE separation!)

% Make a cluster profile and taylor it!
CPR = cluster_profile();
CPR.n_cores_per_node = 4;

% Create a case dispatcher
CPD = case_profile_dispatcher(EDB, CPR);

% Make four simulation profiles
SPR_lb    = simulation_profile();   % Describes lower bounds of parameter values (empty for fixed parameters)
SPR_0     = simulation_profile();   % Describes initial (reference) point (all points filled)
SPR_ub    = simulation_profile();   % Describes upper bounds of parameter values (empty for fixed parameters)
SPR_free  = simulation_profile();   % Highlights free (=1, changed directly by optimizer) and fixed (0, paired to something else)
SPR_index = simulation_profile();   % Describes correspondence to indices in vector passed to optimizer (0 corresponds to copying values from reference case)
% Bernstein coefficients for Cf    correlation
SPR_lb.TA1C  = []   ; SPR_0.TA1C  = 1    ; SPR_ub.TA1C  = []   ; SPR_free.TA1C  = 0    ; SPR_index.TA1C  = 0    ;
SPR_lb.TA2C  = 0.8  ; SPR_0.TA2C  = 1    ; SPR_ub.TA2C  = 1.3  ; SPR_free.TA2C  = 1    ; SPR_index.TA2C  = 1    ;
SPR_lb.TA3C  = 0.8  ; SPR_0.TA3C  = 1    ; SPR_ub.TA3C  = 1.3  ; SPR_free.TA3C  = 1    ; SPR_index.TA3C  = 2    ;
SPR_lb.TA4C  = 0.8  ; SPR_0.TA4C  = 1    ; SPR_ub.TA4C  = 1.3  ; SPR_free.TA4C  = 1    ; SPR_index.TA4C  = 3    ;
SPR_lb.TA5C  = 0.8  ; SPR_0.TA5C  = 1    ; SPR_ub.TA5C  = 1.3  ; SPR_free.TA5C  = 1    ; SPR_index.TA5C  = 4    ;
SPR_lb.TA6C  = []   ; SPR_0.TA6C  = 1    ; SPR_ub.TA6C  = []   ; SPR_free.TA6C  = 0    ; SPR_index.TA6C  = 4    ;
% Bernstein coefficients for Hstar correlation
SPR_lb.TA1H  = []   ; SPR_0.TA1H  = 1    ; SPR_ub.TA1H  = []   ; SPR_free.TA1H  = 0    ; SPR_index.TA1H  = 0    ;
SPR_lb.TA2H  = 0.9  ; SPR_0.TA2H  = 1    ; SPR_ub.TA2H  = 1.1  ; SPR_free.TA2H  = 1    ; SPR_index.TA2H  = 5    ;
SPR_lb.TA3H  = 0.9  ; SPR_0.TA3H  = 1    ; SPR_ub.TA3H  = 1.1  ; SPR_free.TA3H  = 1    ; SPR_index.TA3H  = 6    ;
SPR_lb.TA4H  = 0.9  ; SPR_0.TA4H  = 1    ; SPR_ub.TA4H  = 1.1  ; SPR_free.TA4H  = 1    ; SPR_index.TA4H  = 7    ;
SPR_lb.TA5H  = 0.9  ; SPR_0.TA5H  = 1    ; SPR_ub.TA5H  = 1.1  ; SPR_free.TA5H  = 1    ; SPR_index.TA5H  = 8    ;
SPR_lb.TA6H  = []   ; SPR_0.TA6H  = 1    ; SPR_ub.TA6H  = []   ; SPR_free.TA6H  = 0    ; SPR_index.TA6H  = 8    ;
% Parametrization Properties
SPR_lb.THMIN = []   ; SPR_0.THMIN = 1    ; SPR_ub.THMIN = []   ; SPR_free.THMIN = 0    ; SPR_index.THMIN = 0    ;
SPR_lb.THMAX = []   ; SPR_0.THMAX = 5    ; SPR_ub.THMAX = []   ; SPR_free.THMAX = 0    ; SPR_index.THMAX = 0    ;
SPR_lb.TDCF  = []   ; SPR_0.TDCF  = 0.004; SPR_ub.TDCF  = []   ; SPR_free.TDCF  = 0    ; SPR_index.TDCF  = 0    ;


% Make Global Function Wrapper
GFW = global_function_wrapper(CPD, SPR_lb, SPR_0, SPR_ub, SPR_free, SPR_index);
% Process its initial fields
GFW.process_initial_fields();
% And make intial reference
GFW.make_initial_point_accuracy_profile();

% Now compute a case (we should get scalar_accuracy_metric = 1.9952 with EDB.EC_cell = EC_cell_full(1:10);)
%[scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE, SPR] = accuracy_profiles_from_parameters(GFW, GFW.x_ub);
% Try function handle
scalar_accuracy_metric = scalar_metric_from_parameters(GFW, [GFW.x_0 ; GFW.x_ub]);

% Make function handle for minimization algorithm
min_fun_handle = @(x) scalar_metric_from_parameters(GFW, x);


% Now minimize!
options_grad = optimoptions(@fmincon);
options_grad.FiniteDifferenceStepSize   = 0.0005;
options_grad.FiniteDifferenceType       = 'central';
options_grad.Display                    = 'iter-detailed';

[xgrad,fval_grad,exitflag_grad,output_grad] = fmincon(min_fun_handle ,GFW.x_0,[],[],[],[],GFW.x_lb,GFW.x_ub, [], options_grad);

% Compute population on results
[scalar_accuracy_metric, normalized_accuracy_profile, accuracy_profile, APE, SPR] = accuracy_profiles_from_parameters(GFW, xgrad);

% And plot some stuff
figure(6)
subplot(221)
histogram(GFW.APE_0.cl_accuracy_metric_array); hold on; histogram(accuracy_profile.cl_accuracy_metric_array); legend('Original', 'Optimized'); xlabel('Cl RMS error'); ylabel('# of Airfoils');
title(['C_l accuracy=' num2str(normalized_accuracy_profile.cl_global_accuracy)])
subplot(222)
histogram(GFW.APE_0.cm_accuracy_metric_array); hold on; histogram(accuracy_profile.cm_accuracy_metric_array); legend('Original', 'Optimized'); xlabel('Cm RMS error'); ylabel('# of Airfoils');
title(['C_m accuracy=' num2str(normalized_accuracy_profile.cm_global_accuracy)])
subplot(223)
histogram(GFW.APE_0.cd_accuracy_metric_array, 6); hold on; histogram(accuracy_profile.cd_accuracy_metric_array, 6); legend('Original', 'Optimized'); xlabel('Cd RMS error'); ylabel('# of Airfoils');
title(['C_d accuracy=' num2str(normalized_accuracy_profile.cd_global_accuracy)])
print(['./fig/' , mfilename(),'_', datestr(now(), 30)], '-dpdf') 

%Save results
save([mfilename(),'_' datestr(now(), 30)])

% Now make some more pictures
APE_0 = GFW.APE_0;

for n_foil = 1:length(APE.EC_cell)
    APE.EC_cell{n_foil}.plot()
    figure(1)
    subplot(121); hold on
    plot(  APE.results_list{n_foil}.ap.cd_alpha(  APE.results_list{n_foil}.ap.alpha_range)*1e4,   APE.results_list{n_foil}.ap.cl_alpha(  APE.results_list{n_foil}.ap.alpha_range))
    plot(APE_0.results_list{n_foil}.ap.cd_alpha(APE_0.results_list{n_foil}.ap.alpha_range)*1e4, APE_0.results_list{n_foil}.ap.cl_alpha(APE_0.results_list{n_foil}.ap.alpha_range))
    axis([0 200 -1.5 2])
    title([APE.EC_cell{n_foil}.airfoil_name ' - Re = ' num2str(APE.EC_cell{n_foil}.Re)])
    legend off
    subplot(122); hold on
    plot(  APE.results_list{n_foil}.ap.alpha_range                                            ,   APE.results_list{n_foil}.ap.cl_alpha(  APE.results_list{n_foil}.ap.alpha_range))
    plot(APE_0.results_list{n_foil}.ap.alpha_range                                            , APE_0.results_list{n_foil}.ap.cl_alpha(APE_0.results_list{n_foil}.ap.alpha_range))
    axis([-15 20 -1.5 2])
    title(['\epsilon^{rel}_{Cd} = ' num2str(normalized_accuracy_profile.cd_accuracy_metric_array(n_foil),2) , '   \epsilon^{rel}_{Cl} = ' num2str(normalized_accuracy_profile.cl_accuracy_metric_array(n_foil),2)]);
    legend('TR824 (\alpha-C_l)', 'TR824 (C_l-C_d)', 'LearneRfoil', 'Original Rfoil', 'Location', 'SouthEast')
    print(['./fig/ML_' APE.EC_cell{n_foil}.airfoil_name '_Re', num2str(APE.EC_cell{n_foil}.Re), '_' , datestr(now(), 30)], '-dpdf')
    close 1
end

figure(2)
legend off
title('parantik\epsilon')
legend('run-over-database-alfa15', 'Location', 'South')
print(['./fig/_0_title' , datestr(now(), 30)], '-dpdf')




%system('pdfunite *.pdf out.pdf')
% options = optimoptions(@simulannealbnd, ...
%                      'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping});
% [x,fval,exitFlag,output] = simulannealbnd(min_fun_handle,GFW.x_0,GFW.x_lb,GFW.x_ub,options);

% Now plot some results








