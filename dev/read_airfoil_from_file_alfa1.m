
%% Make a Clean Sheet
clear all; close all; clc

%% Numerical Inputs
N_points              = 160;                              % Number of points to define airfoil shape

%% File Inputs
airfoil_filename      = 'S1223RTL.p160.dat';
%airfoil_filename      = 'RTL50.clean.p160.dat';
%airfoil_filename      = 'MH60_clean.p160.dat';
airfoil_folder        = '/Users/gael/Desktop/eKite/airfoils/';
airfoil_export_prefix = 'S1223RTL';
%airfoil_export_prefix = 'RTL50';
%airfoil_export_prefix = 'MH60';
airfoil_export_suffix = '.air';


%% System Work necessary to start
% Start by adding the source folders to the matlab path
fs = filesep();      % Folder separator is OS dependent
addpath([cd fs 'src']);
addpath([cd fs 'user_src']);
addpath([cd fs 'gui']);

% Create System Context Object and Set Context
SC = system_context; SC.N_cores = 1; SC.set_context;

%% Instanciate geometric manipulations classes
% Parametrization objects for upper and lower side
p_upper = parametrization( 'cst_upper'  , 8);
p_lower = parametrization( 'cst_lower'  , 8);
 
% Fix number of non-CST-shape parameters (1 for TE thickness)
N_dummy_parameters = 1;

% Shape Definition Objects using previously defined parametrizations 
SD_0    = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'cst88');

% Make Shape Fit Object
SF = shape_fit_cst(SD_0, 'cst88');

%% Load Airfoil
fileToRead = [airfoil_folder ,  airfoil_filename];
x_airfoil = SF.get_parameters_from_file(fileToRead);

%% Regenerate Coordinates of Original Airfoil
% First, for unflapped airfoils
[tx_0, tz_0] = SD_0.generate_coordinates(160, x_airfoil);