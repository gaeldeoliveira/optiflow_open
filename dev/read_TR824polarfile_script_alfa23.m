% Script for Reading a Polar File from TR824
% 
%   alfa5 tested manually on 131.ALL and 132.ALL and exhaustively for
%   chapters A and B. Chapter has slightly different filename conventions
%   but we will not use that data explicitly.
%
%   TODO: Check theta sign and power of scale factor!
%   CARE: There is a dangerous warning handler below, just for now!
%#ok<*ASGLU>
%

% Make experimental case database
EDB = experimental_case_database();

% Set Name of Airfoil and Dataset Files
airfoil_folder   = './data/experimental/NACA/TR824-Airfoils/';
dataset_folder   = './data/experimental/NACA/TR824-Digitized/';

% % Now read some data into the system

%% NACA 4 Series (Symmetric)
% Set friendly name of airfoil 
airfoil_name     = 'NACA0006';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca0006.air';
dataset_filename = '131.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA0009';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca0009.air';
dataset_filename = '132.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% NACA 4 Series (14xx)
% Set friendly name of airfoil 
airfoil_name     = 'NACA1408';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca1408.air';
dataset_filename = '133.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA1410';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca1410.air';
dataset_filename = '134.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA1412';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca1412.air';
dataset_filename = '135.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% NACA 4 Series (24xx)
% Set friendly name of airfoil 
airfoil_name     = 'NACA2412';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca2412.air';
dataset_filename = '136.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA2415';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca2415.air';
dataset_filename = '137.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA2418';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca2418.air';
dataset_filename = '138.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA2421';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca2421.air';
dataset_filename = '139.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA2424';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca2424.air';
dataset_filename = '140.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% NACA 4 Series (44xx)
% Set friendly name of airfoil 
airfoil_name     = 'NACA4412';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca4412.air';
dataset_filename = '141.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA4415';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca4415.air';
dataset_filename = '142.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA4418';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca4418.air';
dataset_filename = '143.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA4421';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca4421.air';
dataset_filename = '144.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA4424';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca4424.air';
dataset_filename = '145.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% NACA 4 Series (44xx)
% Set friendly name of airfoil 
airfoil_name     = 'NACA23012';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca23012.air';
dataset_filename = '146.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA23015';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca23015.air';
dataset_filename = '147.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA23018';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca23018.air';
dataset_filename = '148.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA23021';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca23021.air';
dataset_filename = '149.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA23024';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca23024.air';
dataset_filename = '150.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% List of cases (same order and grouping as in present script)
%  131.all  NACA  0006
%  132.all  NACA  0009

%  133.all  NACA  1408
%  134.all  NACA  1410
%  135.all  NACA  1412

%  136.all  NACA  2412
%  137.all  NACA  2415
%  138.all  NACA  2418
%  139.all  NACA  2421
%  140.all  NACA  2424

%  141.all  NACA  4412
%  142.all  NACA  4415
%  143.all  NACA  4418
%  144.all  NACA  4421
%  145.all  NACA  4424

%  146.all  NACA  23012
%  147.all  NACA  23015
%  148.all  NACA  23018
%  149.all  NACA  23021
%  150.all  NACA  23024

