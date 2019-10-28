%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS IS THE "FIRST POST-TORQUE2018 ABSTRACT" Update 
%       Create a new EDB (right now they are working as snapshots)
%       Add NACA 6 series (63,64,65,66) cases
%
% (IMPORTANT CORRECTION: sqrt were not taken in the right place in
%  editions before and up to TORQUE2018 abstract. Corrected 21/02/2018 )
%
%  Narrow down to NACA 4 series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean environment
clear all; close all; clc; lg.clr(mfilename);

% Define script scope with flags
flag_initialization                 = true;
flag_make_new_rng_state             = true;
flag_gradient_descent_deterministic = false;
flag_gradient_descent_stochastic    = false;
flag_postprocess                    = false;

% Define name of state files
initialization_file                 = 'run_over_database_alfa21_initialization_20180223T182432.mat';
rng_state_file                      = 'run_over_database_alfa21_rngstate_20180223T182433.mat';
gradient_file                       = 'run_over_database_alfa21_gradient_20180224T024706.mat';
% Define name of data source files
database_file                       = './data/experimental/EDB_317665fb-cd39-4321-b8cf-124f15d6f060.mat';


% Define operation environment
n_cores_per_node                    = 32;

% Define gradient descent parameters
validation_fraction                 = 0.4;

if flag_initialization == true    
    % System Work necessary to start
    fs = filesep();                 % Folder separator is OS dependent
    addpath([cd fs 'src']);         % Add optimizer source code folder
    addpath([cd fs 'user_src']);    % Add optimizer user sources (cost functions, etc) folder
    addpath([cd fs 'gui']);         % Add gui source code folder
    addpath([cd fs 'dev']);         % Add devsource code folder
    warning('off','MATLAB:DELETE:FileNotFound');
    warning('off','MATLAB:MKDIR:DirectoryExists');
        
    % Get a database
    loaded_datase = load(database_file);
    % Make it present
    EDB = loaded_datase.EDB;
    % Deprecate this (new data only has reasonable cases!)
    % % % Restrict to a subset (for faster prototyping)
    % % EC_cell_full = EDB.EC_cell;
    % % %EDB.EC_cell = EC_cell_full(1:10);
    % % EDB.EC_cell = EC_cell_full(13:60); % Very thin airfoils (first ones, 0006, 0009, 1408 and 1410 lead to issues (convergence?), probably due to LE separation!)
    %
    % Restrict to a subset (for faster prototyping)
    % EC_cell_full = EDB.EC_cell;
    % EDB.EC_cell = EC_cell_full(1:10);
    % Filter out number 135 (NACA 65215 at Re9e6) due to repetition in cd
    % cu
    EC_cell_full = EDB.EC_cell;
    EDB.EC_cell = EDB.EC_cell(2:33);
    
    % Make a cluster profile and taylor it!
    CPR = cluster_profile();
    CPR.n_cores_per_node = n_cores_per_node;
    
    % Create a case dispatcher
    CPD = case_profile_dispatcher(EDB, CPR);
    
    % Make four simulation profiles
    SPR_lb    = simulation_profile();   % Describes lower bounds of parameter values (empty for fixed parameters)
    SPR_0     = simulation_profile();   % Describes initial (reference) point (all points filled)
    SPR_ub    = simulation_profile();   % Describes upper bounds of parameter values (empty for fixed parameters)
    SPR_free  = simulation_profile();   % Highlights free (=1, changed directly by optimizer) and fixed (0, paired to something else)
    SPR_index = simulation_profile();   % Describes correspondence to indices in vector passed to optimizer (0 corresponds to copying values from reference case)
    % Bernstein coefficients for Cf    correlation
    SPR_lb.TA1C  = []   ; SPR_0.TA1C  = 1    ; SPR_ub.TA1C  = []   ; SPR_free.TA1C  = 0    ; SPR_index.TA1C  = 1    ;
    SPR_lb.TA2C  = 0.8  ; SPR_0.TA2C  = 1    ; SPR_ub.TA2C  = 1.3  ; SPR_free.TA2C  = 1    ; SPR_index.TA2C  = 1    ;
    SPR_lb.TA3C  = 0.8  ; SPR_0.TA3C  = 1    ; SPR_ub.TA3C  = 1.3  ; SPR_free.TA3C  = 1    ; SPR_index.TA3C  = 2    ;
    SPR_lb.TA4C  = 0.8  ; SPR_0.TA4C  = 1    ; SPR_ub.TA4C  = 1.3  ; SPR_free.TA4C  = 1    ; SPR_index.TA4C  = 3    ;
    SPR_lb.TA5C  = 0.8  ; SPR_0.TA5C  = 1    ; SPR_ub.TA5C  = 1.3  ; SPR_free.TA5C  = 1    ; SPR_index.TA5C  = 4    ;
    SPR_lb.TA6C  = []   ; SPR_0.TA6C  = 1    ; SPR_ub.TA6C  = []   ; SPR_free.TA6C  = 0    ; SPR_index.TA6C  = 0    ;
    % Bernstein coefficients for Hstar correlation
    SPR_lb.TA1H  = []   ; SPR_0.TA1H  = 1    ; SPR_ub.TA1H  = []   ; SPR_free.TA1H  = 0    ; SPR_index.TA1H  = 5    ;
    SPR_lb.TA2H  = 0.9  ; SPR_0.TA2H  = 1    ; SPR_ub.TA2H  = 1.1  ; SPR_free.TA2H  = 1    ; SPR_index.TA2H  = 5    ;
    SPR_lb.TA3H  = 0.9  ; SPR_0.TA3H  = 1    ; SPR_ub.TA3H  = 1.1  ; SPR_free.TA3H  = 1    ; SPR_index.TA3H  = 6    ;
    SPR_lb.TA4H  = 0.9  ; SPR_0.TA4H  = 1    ; SPR_ub.TA4H  = 1.1  ; SPR_free.TA4H  = 1    ; SPR_index.TA4H  = 7    ;
    SPR_lb.TA5H  = 0.9  ; SPR_0.TA5H  = 1    ; SPR_ub.TA5H  = 1.1  ; SPR_free.TA5H  = 1    ; SPR_index.TA5H  = 8    ;
    SPR_lb.TA6H  = []   ; SPR_0.TA6H  = 1    ; SPR_ub.TA6H  = []   ; SPR_free.TA6H  = 0    ; SPR_index.TA6H  = 0    ;
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
    disp(['Scalar accuracy metric x_0 : ' , num2str(scalar_accuracy_metric(1))]);
    disp(['Scalar accuracy metric x_ub: ' , num2str(scalar_accuracy_metric(2))]);
    
    % Make function handle for minimization algorithm
    min_fun_handle           = @(x) scalar_metric_from_parameters(             GFW, x);
    min_fun_handle_minibatch = @(x) scalar_metric_from_parameters_on_minibatch(GFW, x);
        
    % Make a filename to save current initialization state to file
    initialization_savename = [mfilename() , '_initialization_' , datestr(now(), 30)];
    % Save it to file
    save(initialization_savename);
    % Write a log message
    lg.msg(lg, ['Saved intialization to file : ' , initialization_savename , '.mat'])
else
    % If you didn't compute a new initialization, load it from file!
    load(initialization_file);
    % And report that we loaded it
    lg.msg(lg, ['Loaded initialization from file : ' , initialization_file]);
end

if flag_make_new_rng_state == true
    % Retrieve current state of Random Number Generator
    rngstate = rng();    
    % Make a filename
    rngstate_savename = [mfilename() , '_rngstate_' , datestr(now(), 30)];
    % Save state to file
    save(rngstate_savename, 'rngstate');
    % Write down
    lg.msg(lg, ['Saved RNG state to file : ' , rngstate_savename , '.mat'])
else
    % Load it from file
    loaded_rngstate = load(rng_state_file);
    % Make it accessible
    rngstate = loaded_rngstate.rngstate;
    % And report that we loaded it
    lg.msg(lg, ['Loaded RNG state from file : ' , rng_state_file]);
end
    

if flag_gradient_descent_deterministic == true
%     % Inform we are here, and prepare to split reference dataset
%     lg.msg(lg, 'Enter deteministic gradient module. Split reference dataset');
%     % Define index of reference dataset
     index_reference_dataset = 1:length(EDB.EC_cell); % 1:length(EC_cell_full)
%     
%     % Partition reference dataset with static method of SGD library
%     % (ignore folding: fold_lenght = 1, no monitoring set needed: monitoring_fraction = 0)
%     [index_full_train_subset, ~, index_validation_subset] = ...
%         sgd.shuffle_partition_and_fold_to_length(index_reference_dataset, 1, 0.0, validation_fraction,rngstate);
%     
     % Count datapoints in reference dataset
     n_ref_datapoints   = 0; for n_EC=1:length(index_reference_dataset); n_ref_datapoints   = n_ref_datapoints   + EDB.EC_cell{index_reference_dataset(n_EC)}.count_datapoints(); end
%     % Count datapoints in training dataset
%     n_train_datapoints = 0; for n_EC=1:length(index_full_train_subset); n_train_datapoints = n_train_datapoints + EDB.EC_cell{index_full_train_subset(n_EC)}.count_datapoints(); end
%     % Count datapoints in validation dataset
%     n_val_datapoints   = 0; for n_EC=1:length(index_validation_subset); n_val_datapoints   = n_val_datapoints   + EDB.EC_cell{index_validation_subset(n_EC)}.count_datapoints(); end
    
    % Log datapoint counts
    lg.msg(lg, ['Size of reference  dataset : ' , num2str(n_ref_datapoints  ) , ' datapoints over ' , num2str(length(index_reference_dataset)) ' experimental cases']);
%     lg.msg(lg, ['Size of training   dataset : ' , num2str(n_train_datapoints) , ' datapoints over ' , num2str(length(index_full_train_subset)) ' experimental cases']);
%     lg.msg(lg, ['Size of validation dataset : ' , num2str(n_val_datapoints  ) , ' datapoints over ' , num2str(length(index_validation_subset)) ' experimental cases']);
%     
%     lg.msg(lg, 'Test innacuracy norm on training and validation set.');
%     min_fun_handle_full_train_set = @(x) scalar_metric_from_parameters_on_minibatch(GFW, x, index_full_train_subset);
%     min_fun_handle_validation_set = @(x) scalar_metric_from_parameters_on_minibatch(GFW, x, index_validation_subset);
    
    % Define 'classic' graident minimization problem!
    options_grad = optimoptions(@fmincon);
    options_grad.FiniteDifferenceStepSize   = 0.0001;           % Original runs with 0.0005. Screws up at 0.00001 for first coordinate, still works at 0.00005 with EPSV=0=1e-4 (works at least down to 0.00001 with EPSV=1E-6)
    options_grad.FiniteDifferenceType       = 'central';
    options_grad.Display                    = 'iter-detailed';
    % Now minimize!
    lg.msg(lg, 'Start solution of minimization problem.');
    [xgrad,fval_grad,exitflag_grad,output_grad] = fmincon(min_fun_handle ,GFW.x_0,[],[],[],[],GFW.x_lb,GFW.x_ub, [], options_grad);
    lg.msg(lg, 'Done with minimization problem.');
    
    % Compute dataset innacuracy measures with new optimized closures
    [scalar_accuracy_metric           , normalized_accuracy_profile           , accuracy_profile           , APE           , SPR           ] = GFW.accuracy_profiles_from_parameters(xgrad);    
    %[scalar_accuracy_metric_full_train, normalized_accuracy_profile_full_train, accuracy_profile_full_train, APE_full_train, SPR_full_train] = GFW.accuracy_profiles_from_parameters_on_minibatch(xgrad, index_full_train_subset);
    %[scalar_accuracy_metric_validation, normalized_accuracy_profile_validation, accuracy_profile_validation, APE_validation, SPR_validation] = GFW.accuracy_profiles_from_parameters_on_minibatch(xgrad, index_validation_subset);
    
    % Make a filename to save current results (full-state) to file
    gradient_savename = [mfilename() , '_gradient_' , datestr(now(), 30)];
    % Save it to file
    save(gradient_savename);
    % Write a log message
    lg.msg(lg, ['Saved after deterministic gradient descent to file : ' , gradient_savename , '.mat'])
    
end

if flag_postprocess == true
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
        savefig(['./fig/ML_' APE.EC_cell{n_foil}.airfoil_name '_Re', num2str(APE.EC_cell{n_foil}.Re), '_' , datestr(now(), 30)]);
        close 1
    end
    
    figure(2)
    legend off
    title('parantik\epsilon')
    legend(mfilename(), 'Location', 'South')
    print(  ['./fig/' , mfilename() , '_0_title' , datestr(now(), 30)], '-dpdf')
    savefig(['./fig/' , mfilename() , '_0_title' , datestr(now(), 30)]);
    
    %system('pdfunite *.pdf out.pdf')
    % options = optimoptions(@simulannealbnd, ...
    %                      'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping});
    % [x,fval,exitFlag,output] = simulannealbnd(min_fun_handle,GFW.x_0,GFW.x_lb,GFW.x_ub,options);
    
end








