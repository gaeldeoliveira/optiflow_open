% Script for Reading a Polar File from TR824
% 
%   alfa5 tested manually on 131.ALL and 132.ALL and exhaustively for
%   chapters A and B. Chapter has slightly different filename conventions
%   but we will not use that data explicitly.
%
%   TODO: Check theta sign and power of scale factor!

% Make experimental case database
EDB = experimental_case_database();

% Set Name of Airfoil and Dataset Files
airfoil_folder   = './data/experimental/NACA/TR824-Airfoils/';
dataset_folder   = './data/experimental/NACA/TR824-Digitized/';

% % Now read some data into the system

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






% Plotting of processed polars
figure(1)
processed_polars = EDB.EC_cell{3}.processed_polars;
subplot(122)
plot(processed_polars.al_alcl_polar   , processed_polars.cl_alcl_polar   , 'x-'); hold on;
plot(processed_polars.al_alclcd_polar , processed_polars.cl_alclcd_polar , 'o-'); grid on;
xlabel('\alpha (deg)'); ylabel('C_l');
subplot(121)
plot(processed_polars.cd_clcd_polar   , processed_polars.cl_clcd_polar   , 'x-'); hold on;
plot(processed_polars.cd_alclcd_polar , processed_polars.cl_alclcd_polar , 'o-'); grid on;
xlabel('C_d')         ; ylabel('C_l');
legend('Sorted Data' , 'Reinterpolated Data', 'Location', 'East');

figure(2)
plot(airfoil_description.tx_coordinates    , airfoil_description.tz_coordinates          ); hold on; 
plot(airfoil_description.tx_coordinates_raw, airfoil_description.tz_coordinates_raw, '--.'); grid on; 
[tx, tz] = EDB.SD.generate_coordinates(200, airfoil_description.x);
plot(tx, tz, '-.'); grid on; 
xlabel('x/c'); ylabel('y/c'); axis equal;
legend('Standardized Airfoil Shape' , 'Raw Airfoil Shape', 'Regenerated Airfoil Shape');
print([airfoil_filename '-polar'], '-dpdf')




% Create interpolants 







% 
%Re_dataset_01 = 3e6;
%cl_polar_alpha_data_dataset_01

%dataset_01_alpha_cl = ;
%dataset_01_alpha_for_cl = ;
%dataset_01_alpha_for_cl = ;

