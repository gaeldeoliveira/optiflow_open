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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS IS THE "FIRST POST-TORQUE2018 ABSTRACT" Update 
%       Create a new EDB (right now they are working as snapshots)
%       Add NACA 6 series (63,64,65,66) cases
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Make experimental case database
EDB = experimental_case_database();

% Set Name of Airfoil and Dataset Files
airfoil_folder   = './data/experimental/NACA/TR824-Airfoils/';
dataset_folder   = './data/experimental/NACA/TR824-Digitized/';

% % Now read some data into the system

%% NACA 4 Series (Symmetric, don't take airfoils below 12% thickness to avoid LE separation issues)
% % Set friendly name of airfoil 
% airfoil_name     = 'NACA0006';
% % Set filenames for airfoil and dataset Files
% airfoil_filename = 'naca0006.air';
% dataset_filename = '131.ALL';
% % Load Make and Add airfoil Description
% [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);
% 
% % Set friendly name of airfoil 
% airfoil_name     = 'NACA0009';
% % Set filenames for airfoil and dataset Files
% airfoil_filename = 'naca0009.air';
% dataset_filename = '132.ALL';
% % Load Make and Add airfoil Description
% [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);
% 
% %% NACA 4 Series (14xx)
% % Set friendly name of airfoil 
% airfoil_name     = 'NACA1408';
% % Set filenames for airfoil and dataset Files
% airfoil_filename = 'naca1408.air';
% dataset_filename = '133.ALL';
% % Load Make and Add airfoil Description
% [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);
% 
% % Set friendly name of airfoil 
% airfoil_name     = 'NACA1410';
% % Set filenames for airfoil and dataset Files
% airfoil_filename = 'naca1410.air';
% dataset_filename = '134.ALL';
% % Load Make and Add airfoil Description
% [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% NACA 4 Series (14xx)

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

%% Now add the 6 series airfoils

%% 63 Series 12 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA63012';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-012.gnu';
dataset_filename = '163.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA63212';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-212.gnu';
dataset_filename = '164.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA63412';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-412.gnu';
dataset_filename = '165.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);


%% 63 Series 15 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA63015';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-015.gnu';
dataset_filename = '166.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA63215';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-215.gnu';
dataset_filename = '167.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA63415';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-415.gnu';
dataset_filename = '168.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA63615';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-615.gnu';
dataset_filename = '169.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);


%% 63 Series 18 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA63018';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-018.gnu';
dataset_filename = '170.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA63218';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-218.gnu';
dataset_filename = '171.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA63418';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-418.gnu';
dataset_filename = '172.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA63618';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'NACA-63618.lc';
dataset_filename = '173.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 63 Series 21 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA63021';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-021.gnu';
dataset_filename = '174.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA63221';
% Set filenames for airfoil and dataset Files
airfoil_filename = '63-221.gnu';
dataset_filename = '175.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA63421';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'NACA-63421.lc';
dataset_filename = '176.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% Now add the 63 series airfoils

%  163.all  NACA  63(sub)1-012  OK
%  164.all  NACA  63(sub)1-212  OK
%  165.all  NACA  63(sub)1-412  OK

%  166.all  NACA  63(sub)2-015  OK
%  167.all  NACA  63(sub)2-215  OK
%  168.all  NACA  63(sub)2-415  OK
%  169.all  NACA  63(sub)2-615  OK

%  170.all  NACA  63(sub)3-018  OK
%  171.all  NACA  63(sub)3-218  OK
%  172.all  NACA  63(sub)3-418  OK
%  173.all  NACA  63(sub)3-618  OK

%  174.all  NACA  63(sub)4-021  OK
%  175.all  NACA  63(sub)4-221  OK
%  176.all  NACA  63(sub)4-421  OK

%% 64 Series 12 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA64012';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-012.gnu';
dataset_filename = '185.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA64112';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-112.gnu';
dataset_filename = '186.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA64212';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-212.gnu';
dataset_filename = '187.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA64412';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-412.gnu';
dataset_filename = '188.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 64 Series 15 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA64015';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-015.gnu';
dataset_filename = '189.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA64215';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-215.gnu';
dataset_filename = '190.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA64415';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-415.gnu';
dataset_filename = '191.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 64 Series 18 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA64018';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-018.gnu';
dataset_filename = '192.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA64218';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-218.gnu';
dataset_filename = '193.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA64418';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-418.gnu';
dataset_filename = '194.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA64618';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'NACA-64618.lc';
dataset_filename = '195.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 64 Series 21 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA64021';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-021.gnu';
dataset_filename = '196.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA64221';
% Set filenames for airfoil and dataset Files
airfoil_filename = '64-221.gnu';
dataset_filename = '197.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% % Set friendly name of airfoil 
airfoil_name     = 'NACA64421';
% Set filenames for airfoil and dataset Files
airfoil_filename = 'naca634421_coord.txt';
dataset_filename = '198.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% Now add the 64 series airfoils
%  
%  185.all  NACA  64(sub)1-012  OK
%  186.all  NACA  64(sub)1-112  OK
%  187.all  NACA  64(sub)1-212  OK
%  188.all  NACA  64(sub)1-412  OK

%  189.all  NACA  64(sub)2-015  OK
%  190.all  NACA  64(sub)2-215  OK
%  191.all  NACA  64(sub)2-415  OK

%  192.all  NACA  64(sub)3-018  OK
%  193.all  NACA  64(sub)3-218  OK
%  194.all  NACA  64(sub)3-418  OK
%  195.all  NACA  64(sub)3-618  OK

%  196.all  NACA  64(sub)4-021  OK
%  197.all  NACA  64(sub)4-221  OK
%  198.all  NACA  64(sub)4-421  OK

%% 65 Series 12 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA65012';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-012.gnu';
dataset_filename = '210.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA65212';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-212.gnu';
dataset_filename = '211.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA65412';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-412.gnu';
dataset_filename = '214.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 65 Series 15 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA65015';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-015.gnu';
dataset_filename = '215.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA65215';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-215.gnu';
dataset_filename = '216.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA65415';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-415.gnu';
dataset_filename = '217.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 65 Series 18 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA65018';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-018.gnu';
dataset_filename = '219.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA65218';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-218.gnu';
dataset_filename = '222.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA65418';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-418.gnu';
dataset_filename = '223.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% % Set friendly name of airfoil 
% airfoil_name     = 'NACA65618';
% % Set filenames for airfoil and dataset Files
% airfoil_filename = '65-618.gnu';
% dataset_filename = '225.ALL';
% % Load Make and Add airfoil Description
% [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);


%% 65 Series 21 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA65021';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-021.gnu';
dataset_filename = '227.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA65221';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-221.gnu';
dataset_filename = '228.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA65421';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-421.gnu';
dataset_filename = '229.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% Now add the 65 series airfoils
 
%  210.all  NACA  65(sub)1-012  OK
%  211.all  NACA  65(sub)1-212  OK
%  214.all  NACA  65(sub)1-412  OK

%  215.all  NACA  65(sub)2-015  OK
%  216.all  NACA  65(sub)2-215  OK
%  217.all  NACA  65(sub)2-415  OK
% 
%  219.all  NACA  65(sub)3-018  OK
%  222.all  NACA  65(sub)3-218  OK
%  223.all  NACA  65(sub)3-418  OK
%  225.all  NACA  65(sub)3-618  OK
% 
%  227.all  NACA  65(sub)4-021  OK
%  228.all  NACA  65(sub)4-221  OK
%  229.all  NACA  65(sub)4-421  OK

%% 66 Series 12 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA66012';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-012.gnu';
dataset_filename = '252.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA66212';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-212.gnu';
dataset_filename = '253.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 66 Series 15 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA66015';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-015.gnu';
dataset_filename = '254.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA66215';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-215.gnu';
dataset_filename = '255.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA66415';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-415.gnu';
dataset_filename = '256.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 66 Series 18 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA66018';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-018.gnu';
dataset_filename = '257.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA66218';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-218.gnu';
dataset_filename = '258.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil 
airfoil_name     = 'NACA66418';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-418.gnu';
dataset_filename = '259.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% 66 Series 21 Percent Thick Airfoils
% Set friendly name of airfoil 
airfoil_name     = 'NACA66021';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-021.gnu';
dataset_filename = '260.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (point missmatch near stall between two RE curves?)
airfoil_name     = 'NACA66221';
% Set filenames for airfoil and dataset Files
airfoil_filename = '66-221.gnu';
dataset_filename = '261.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% Now add the 66 series airfoils
% 
%  252.all  NACA  66(sub)1-012  OK
%  253.all  NACA  66(sub)1-212  OK

%  254.all  NACA  66(sub)2-015  OK
%  255.all  NACA  66(sub)2-215  OK
%  256.all  NACA  66(sub)2-415  OK

%  257.all  NACA  66(sub)3-018  OK
%  258.all  NACA  66(sub)3-218  OK
%  259.all  NACA  66(sub)3-418  OK

%  260.all  NACA  66(sub)4-021  OK
%  261.all  NACA  66(sub)4-221  OK

%% 66 Series Airfoils with Custom Camberline 

% Set friendly name of airfoil 
airfoil_name     = 'NACA65212 a=0.6';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-212a06.gnu';
dataset_filename = '213.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (soft stall, sensitive to reynolds, very nice!)
airfoil_name     = 'NACA65415 a=0.5';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-415a05.gnu';
dataset_filename = '218.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% Set friendly name of airfoil (soft stall, sensitive to reynolds, very nice!)
airfoil_name     = 'NACA65418 a=0.5';
% Set filenames for airfoil and dataset Files
airfoil_filename = '65-418a05.gnu';
dataset_filename = '224.ALL';
% Load Make and Add airfoil Description
[EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% % Set friendly name of airfoil (soft stall, sensitive to reynolds, very nice!)
% airfoil_name     = 'NACA65618 a=0.5';
% % Set filenames for airfoil and dataset Files
% airfoil_filename = '65-618a05.gnu';
% dataset_filename = '226.ALL';
% % Load Make and Add airfoil Description
% [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

% % Set friendly name of airfoil (soft stall, sensitive to reynolds, very nice!)
% airfoil_name     = 'NACA65421 a=0.5';
% % Set filenames for airfoil and dataset Files
% airfoil_filename = '65-421a05.gnu';
% dataset_filename = '230.ALL';
% % Load Make and Add airfoil Description
% [EC_re3e6clean , EC_re6e6clean , EC_re9e6clean, airfoil_description] = TR824_reader.load_make_and_add_experimental_cases(EDB, airfoil_name, [airfoil_folder , airfoil_filename], [dataset_folder, dataset_filename]);

%% Now add the custom camberline 65 series airfoils
%  
%  213.all  NACA  65(sub)1-212, a = 0.6     OK
%  218.all  NACA  65(sub)2-415, a = 0.5     OK
%  224.all  NACA  65(sub)3-418, a = 0.5     OK
%  226.all  NACA  65(sub)3-618, a = 0.5     Keep Out: Coordinate Fitting Issue
%  230.all  NACA  65(sub)4-421, a = 0.5     Keep Out: Coordinate Fitting Issue

%% Other rulled out airfoils
%   NACA63618	FAIL		|	FAIL
%   NACA63421	FAIL		|	FAIL
%   NACA64618	FAIL		|	FAIL
%   NACA64421	FAIL		|	FAIL
%   NACA65618	FAIL		|	FAIL

%% Verify fit
plot(EDB.EC_cell{end}.al_alcl_polar_std, EDB.EC_cell{end}.cl_alcl_polar_std); grid on; hold on; plot(EDB.EC_cell{end-1}.al_alcl_polar_std, EDB.EC_cell{end-1}.cl_alcl_polar_std); plot(EDB.EC_cell{end-2}.al_alcl_polar_std, EDB.EC_cell{end-2}.cl_alcl_polar_std); grid on; hold on

for n_EC_cell=1:3:length(EDB.EC_cell)
    EDB.SD.plot_airfoil(EDB.EC_cell{n_EC_cell}.airfoil_description.x); hold on; 
    plot(EDB.EC_cell{n_EC_cell}.airfoil_description.tx_coordinates_raw, EDB.EC_cell{n_EC_cell}.airfoil_description.tz_coordinates_raw, '--');
    plot(EDB.EC_cell{n_EC_cell}.airfoil_description.tx_coordinates, EDB.EC_cell{n_EC_cell}.airfoil_description.tz_coordinates, '--'); hold off; 
    title([EDB.EC_cell{n_EC_cell}.airfoil_name, ' - reconstruction with ', 'CST' , num2str(EDB.SD.parametrization_handles.upper.order) , '' , num2str(EDB.SD.parametrization_handles.lower.order)])
    print('-dpdf', ['foil_reconstruction_', EDB.EC_cell{n_EC_cell}.airfoil_name, '.pdf'])
    axis([0 0.12 -.06 .06]); %axis([0 0.2 -.1 .1])
    print('-dpdf', ['foil_reconstruction_LEzoom_', EDB.EC_cell{n_EC_cell}.airfoil_name, '.pdf'])
    close all
end

for n_EC_cell=1:3:length(EDB.EC_cell)
    disp(EDB.EC_cell{n_EC_cell}.airfoil_name)
end

%% Now save data 
%save(FILENAME,VARIABLES)
db_filename = ['EDB_' , EDB.uuid , '.mat'];
db_folder   = ['./data/experimental/'    ];

save([db_folder , db_filename],'EDB')

%% Summary of Fitting Quality
% NACA1412	OK          |	OK
% NACA2412	OK          |	OK
% NACA2415	OK          |	OK
% NACA2418	OKish		|	OK
% NACA2421	OKish		|	OK
% NACA2424	OKish		|	OK
% NACA4412	OK          |	OK
% NACA4415	OKish		|	OK
% NACA4418	OK          |	OK
% NACA4421	OKish		|	OK
% NACA4424	OK          |	OK
% NACA23012	OK          |	OK
% NACA23015	OK          |	OK
% NACA23018	OK          |	OK
% NACA23021	OK          |	OK	
% NACA23024	OK          |	OK
% NACA63012	OK          |	OK
% NACA63212	OKish		|	OK
% NACA63412	OKish		|	OK
% NACA63015	OK          |	OK
% NACA63215	OKish		|	OK
% NACA63415	OKishish	|	OK
% NACA63615	OKishish	|	OK
% NACA63018	OK          |	OK
% NACA63218	OKishish	|	OK
% NACA63418	OKishish	|	OK
% NACA63618	FAIL		|	FAIL
% NACA63021	OK          |	OK
% NACA63221	OKishsishish|	OK
% NACA63421	FAIL		|	FAIL
% NACA64012	OK          |	OK
% NACA64112	OK          |	OK
% NACA64212	OK          |	OK
% NACA64412	OKish		|	OK
% NACA64015	OK          |	OK
% NACA64215	OK          |	OK
% NACA64415	OKish		|	OK
% NACA64018	OK          |	OK
% NACA64218	OKish		|	OK
% NACA64418	OKishish	|	OK
% NACA64618	FAIL		|	FAIL
% NACA64021	OK          |	OK
% NACA64221	OKish		|	OK
% NACA64421	FAIL		|	FAIL
% NACA65012	OK          |	OK
% NACA65212	OK          |	OK
% NACA65412	OKish		|	OK
% NACA65015	OK          |	OK
% NACA65215	OK          |	OK
% NACA65415	OKish		|	OK
% NACA65018	OK          |	OK
% NACA65218	OK          |	OK
% NACA65418	OKishishish	|	OK
% NACA65618	FAIL		|	FAIL
% NACA65021	OK          |	OK
% NACA65221	OK          |	OK
% NACA65421	OKishishish	|	OK
% NACA66012	OK          |	OK
% NACA66212	OK          |	OK
% NACA66015	OK          |	OK
% NACA66215	OK          |	OK
% NACA66415	OKishish	|	
% NACA66018	OK          |	OK
% NACA66218	OKish		|	OK
% NACA66418	OKishish	|
% NACA66021	OK          |	OK
% NACA66221	OKish		|	OK
% NACA65212 a=0.6	OK	|	OK
% NACA65415 a=0.5	FAIL|	OKish
% NACA65418 a=0.5	FAIL|	OKish

