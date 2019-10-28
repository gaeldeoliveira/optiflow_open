
% Creates a database and loads experimental data into it

% Start from a clean sheet
clear all; close all; clc;

% Make experimental case database
EDB = experimental_case_database();

% Load TR824 airfoil cases
read_TR824polarfile_script_alfa24();

% Load BL Run cases
load_BLdata_cases_alfa2();

% Save database
save(['./data/experimental/EDB' EDB.uuid], 'EDB');
