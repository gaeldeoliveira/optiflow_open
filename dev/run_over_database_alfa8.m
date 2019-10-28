% Clean environment
clear all; close all; clc;

% System Work necessary to start
fs = filesep();                 % Folder separator is OS dependent
addpath([cd fs 'src']);         % Add optimizer source code folder
addpath([cd fs 'user_src']);    % Add optimizer user sources (cost functions, etc) folder
addpath([cd fs 'gui']);         % Add gui source code folder
addpath([cd fs 'dev']);         % Add devsource code folder
warning('off','MATLAB:DELETE:FileNotFound');

% Get a database
loaded_datase = load('./data/experimental/EDB18dfdd03-e00b-4717-8af4-6ba3892db2e5.mat');
% Make it present
EDB = loaded_datase.EDB;

% Make a cluster profile
CPR = cluster_profile();

% Create a case dispatcher
CPD = case_profile_dispatcher(EDB, CPR);

% Make four simulation profiles
SPR_lb    = simulation_profile();   % Describes lower bounds of parameter values (empty for fixed parameters)
SPR_0     = simulation_profile();   % Describes initial (reference) point (all points filled)
SPR_ub    = simulation_profile();   % Describes upper bounds of parameter values (empty for fixed parameters)
SPR_free  = simulation_profile();   % Highlights free (=1, changed directly by optimizer) and fixed (0, paired to something else)
SPR_index = simulation_profile();   % Describes correspondence to indices in vector passed to optimizer (0 corresponds to copying values from reference case)
% Bernstein coefficients for Cf    correlation
SPR_lb.TA1C  = 0.8  ; SPR_0.TA1C  = 1    ; SPR_ub.TA1C  = 1.3  ; SPR_free.TA1C  = 1    ; SPR_index.TA1C  = 1    ;
SPR_lb.TA2C  = 0.8  ; SPR_0.TA2C  = 1    ; SPR_ub.TA2C  = 1.3  ; SPR_free.TA2C  = 1    ; SPR_index.TA2C  = 2    ;
SPR_lb.TA3C  = 0.8  ; SPR_0.TA3C  = 1    ; SPR_ub.TA3C  = 1.3  ; SPR_free.TA3C  = 1    ; SPR_index.TA3C  = 3    ;
SPR_lb.TA4C  = 0.8  ; SPR_0.TA4C  = 1    ; SPR_ub.TA4C  = 1.3  ; SPR_free.TA4C  = 1    ; SPR_index.TA4C  = 4    ;
SPR_lb.TA5C  = []   ; SPR_0.TA5C  = 1    ; SPR_ub.TA5C  = []   ; SPR_free.TA5C  = 0    ; SPR_index.TA5C  = 4    ;
SPR_lb.TA6C  = []   ; SPR_0.TA6C  = 1    ; SPR_ub.TA6C  = []   ; SPR_free.TA6C  = 0    ; SPR_index.TA6C  = 4    ;
% Bernstein coefficients for Hstar correlation
SPR_lb.TA1H  = []   ; SPR_0.TA1H  = 1    ; SPR_ub.TA1H  = []   ; SPR_free.TA1H  = 0    ; SPR_index.TA1H  = 0    ;
SPR_lb.TA2H  = 0.9  ; SPR_0.TA2H  = 1    ; SPR_ub.TA2H  = 1.1  ; SPR_free.TA2H  = 1    ; SPR_index.TA2H  = 5    ;
SPR_lb.TA3H  = 0.9  ; SPR_0.TA3H  = 1    ; SPR_ub.TA3H  = 1.1  ; SPR_free.TA3H  = 1    ; SPR_index.TA3H  = 6    ;
SPR_lb.TA4H  = 0.9  ; SPR_0.TA4H  = 1    ; SPR_ub.TA4H  = 1.1  ; SPR_free.TA4H  = 1    ; SPR_index.TA4H  = 7    ;
SPR_lb.TA5H  = []   ; SPR_0.TA5H  = 1    ; SPR_ub.TA5H  = []   ; SPR_free.TA5H  = 0    ; SPR_index.TA5H  = 7    ;
SPR_lb.TA6H  = []   ; SPR_0.TA6H  = 1    ; SPR_ub.TA6H  = []   ; SPR_free.TA6H  = 0    ; SPR_index.TA6H  = 7    ;
% Parametrization Properties
SPR_lb.THMIN = []   ; SPR_0.THMIN = 1    ; SPR_ub.THMIN = []   ; SPR_free.THMIN = 0    ; SPR_index.THMIN = 0    ;
SPR_lb.THMAX = []   ; SPR_0.THMAX = 5    ; SPR_ub.THMAX = []   ; SPR_free.THMAX = 0    ; SPR_index.THMAX = 0    ;
SPR_lb.TDCF  = []   ; SPR_0.TDCF  = 0.004; SPR_ub.TDCF  = []   ; SPR_free.TDCF  = 0    ; SPR_index.TDCF  = 0    ;


% Make Global Function Wrapper
GFW = global_function_wrapper(CPD, SPR_lb, SPR_0, SPR_ub, SPR_free, SPR_index);
% And process its initial fields
GFW.process_initial_fields();

% Now define a a parameter vector
x = [1.2 1.3 1.4 1.5 0.9 0.8 0.7];

% % And reconstruct an SPR from parameter vector (receive x, return SPR)
% Create virgin simulation profile
SPR = simulation_profile();
% Now fill it in!
for n_field = 1:length(GFW.SPR_fields)
    % Get parameter index of current field
    i_parameter       = GFW.SPR_index.(GFW.SPR_fields{n_field});
    
    % If field is defined from parameter vector (may be non-free, but duplicated)
    if not(i_parameter == 0)
        % Fill it in from parameter
        SPR.(GFW.SPR_fields{n_field}) = x(i_parameter);
    else
        % Duplicate it from template (initial point)
        SPR.(GFW.SPR_fields{n_field}) = GFW.SPR_0.(GFW.SPR_fields{n_field});
    end
end










% Restrict to a subset (for faster prototyping)
EC_cell_full = CPD.EDB.EC_cell;
CPD.EDB.EC_cell = EC_cell_full(1:10);
% Request the evaluation of a simulation profile
[accuracy_profile_0, APE_0] = CPD.dispatch_simulation_profile(SPR_0);
[accuracy_profile_0, APE_lb] = CPD.dispatch_simulation_profile(SPR_lb);
[accuracy_profile_0, APE_ub] = CPD.dispatch_simulation_profile(SPR_ub);

% % Case function part
%EC_cell_airfoil_cases

%[cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper(CPR, SPR, EC_cell_airfoil_cases)




hist(accuracy_profile.cl_accuracy_metric_array, 3)











