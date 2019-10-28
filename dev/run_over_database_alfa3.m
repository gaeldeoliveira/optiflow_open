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

% Make a dummy simulation profile
SPR = [];

% Request the evaluation of a simulation profile
CPD.dispatch_simulation_profile(SPR);

% % Case function part
%EC_cell_airfoil_cases

%[cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper(CPR, SPR, EC_cell_airfoil_cases)
















