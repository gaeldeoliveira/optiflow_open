% Script for Reading a Polar File from TR824
% 
%   alfa5 tested manually on 131.ALL and 132.ALL and exhaustively for
%   chapters A and B. Chapter has slightly different filename conventions
%   but we will not use that data explicitly.


% % Now read some data into the system

% Set Name of File to Read Datasets
dataset_filename = './data/experimental/NACA/TR824-Digitized/132.ALL';

% Load datasets
[datasets] = TR824_reader.read_datasets_from_file(filename);

% % % Process Polars from Datasets
% % Processs polars for Re=3e6 datasets
% Let us now extract a tripplet of datasets
alcl_polar_dataset_re3e6clean     = datasets.A01_data;
alcm_polar_dataset_re3e6clean     = datasets.A07_data;
clcd_polar_dataset_re3e6clean     = datasets.B01_data;
% And get some processed polars out of it
[processed_polars_re3e6clean]     = TR824_reader.process_polar_triplet(alcl_polar_dataset_re3e6clean, alcm_polar_dataset_re3e6clean, clcd_polar_dataset_re3e6clean);

% % Processs polars for Re=6e6 datasets
% Let us now extract a tripplet of datasets
alcl_polar_dataset_re6e6clean     = datasets.A02_data;
alcm_polar_dataset_re6e6clean     = datasets.A08_data;
clcd_polar_dataset_re6e6clean     = datasets.B02_data;
% And get some processed polars out of it
[processed_polars_re6e6clean]     = TR824_reader.process_polar_triplet(alcl_polar_dataset_re6e6clean, alcm_polar_dataset_re6e6clean, clcd_polar_dataset_re6e6clean);

% % Process polars for Re=9e6 datasets
% Let us now extract a tripplet of datasets
alcl_polar_dataset_re9e6clean     = datasets.A03_data;
alcm_polar_dataset_re9e6clean     = datasets.A09_data;
clcd_polar_dataset_re9e6clean     = datasets.B03_data;
% And get some processed polars out of it
[processed_polars_re9e6clean]     = TR824_reader.process_polar_triplet(alcl_polar_dataset_re9e6clean, alcm_polar_dataset_re9e6clean, clcd_polar_dataset_re9e6clean);

% % % Make experimental polar conditions
% % Make experimental polar conditions for Re=3e6 clean polars
% Create structure (object) from template
polar_conditions_re3e6clean = experimental_polar_conditions();
% Fill it in
polar_conditions_re3e6clean.M        = 0;     % Mach number
polar_conditions_re3e6clean.Re       = 3e6;   % Reynolds Number
polar_conditions_re3e6clean.N_crit   = 9;     % Critical Amplification Factor of Tollmien-Schlichted waves (Van Ingen e^N method)
polar_conditions_re3e6clean.xtr_top  = 0.99;  % x/c of forced transition point for top side (1 = free)
polar_conditions_re3e6clean.xtr_bot  = 0.99;  % x/c of forced transition point for bottom side (1 = free)

% % Make experimental polar conditions for Re=6e6 clean polars
% Create structure (object) from template
polar_conditions_re6e6clean = experimental_polar_conditions();
% Fill it in
polar_conditions_re6e6clean.M        = 0;     % Mach number
polar_conditions_re6e6clean.Re       = 6e6;   % Reynolds Number
polar_conditions_re6e6clean.N_crit   = 9;     % Critical Amplification Factor of Tollmien-Schlichted waves (Van Ingen e^N method)
polar_conditions_re6e6clean.xtr_top  = 0.99;  % x/c of forced transition point for top side (1 = free)
polar_conditions_re6e6clean.xtr_bot  = 0.99;  % x/c of forced transition point for bottom side (1 = free)

% % Make experimental polar conditions for Re=9e6 clean polars
% Create structure (object) from template
polar_conditions_re9e6clean = experimental_polar_conditions();
% Fill it in
polar_conditions_re9e6clean.M        = 0;     % Mach number
polar_conditions_re9e6clean.Re       = 9e6;   % Reynolds Number
polar_conditions_re9e6clean.N_crit   = 9;     % Critical Amplification Factor of Tollmien-Schlichted waves (Van Ingen e^N method)
polar_conditions_re9e6clean.xtr_top  = 0.99;  % x/c of forced transition point for top side (1 = free)
polar_conditions_re9e6clean.xtr_bot  = 0.99;  % x/c of forced transition point for bottom side (1 = free)

% % % Nake experimental cases







% Plotting of processed polars
figure(2)
processed_polars = processed_polars_re9e6clean;
subplot(122)
plot(processed_polars.al_alcl_polar   , processed_polars.cl_alcl_polar   , 'x-'); hold on;
plot(processed_polars.al_alclcd_polar , processed_polars.cl_alclcd_polar , 'o-'); grid on;
xlabel('\alpha (deg)'); ylabel('C_l');
subplot(121)
plot(processed_polars.cd_clcd_polar   , processed_polars.cl_clcd_polar   , 'x-'); hold on;
plot(processed_polars.cd_alclcd_polar , processed_polars.cl_alclcd_polar , 'o-'); grid on;
xlabel('C_d')         ; ylabel('C_l');
legend('Sorted Data' , 'Reinterpolated Data', 'Location', 'East')


% Create interpolants 







% 
%Re_dataset_01 = 3e6;
%cl_polar_alpha_data_dataset_01

%dataset_01_alpha_cl = ;
%dataset_01_alpha_for_cl = ;
%dataset_01_alpha_for_cl = ;

