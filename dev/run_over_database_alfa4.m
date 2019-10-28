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

% Make an initial simulation profile
SPR = simulation_profile();
% Bernstein coefficients for Cf    correlation
SPR.TA1C  = 1; SPR.TA2C  = 1; SPR.TA3C  = 1;
SPR.TA4C  = 1; SPR.TA5C  = 1; SPR.TA6C  = 1;
% Bernstein coefficients for Hstar correlation
SPR.TA1H  = 1; SPR.TA2H  = 1; SPR.TA3H  = 1;
SPR.TA4H  = 1; SPR.TA5H  = 1; SPR.TA6H  = 1;
% Parametrization Properties
SPR.THMIN = 1; SPR.THMAX = 5; SPR.TDCF  = 0.004;

% Restrict to a subset
EC_cell_full = CPD.EDB.EC_cell;
CPD.EDB.EC_cell = EC_cell_full(1);
% Request the evaluation of a simulation profile
[accuracy_profile, APE] = CPD.dispatch_simulation_profile(SPR);

% % Case function part
%EC_cell_airfoil_cases

%[cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper(CPR, SPR, EC_cell_airfoil_cases)
















