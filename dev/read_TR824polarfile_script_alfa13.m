% Script for Reading a Polar File from TR824
% 
%   alfa5 tested manually on 131.ALL and 132.ALL and exhaustively for
%   chapters A and B. Chapter has slightly different filename conventions
%   but we will not use that data explicitly.


% % Now read some data into the system

% Set Name of File to Read
filename = './data/experimental/NACA/TR824-Digitized/132.ALL';

[datasets] = TR824_reader.read_datasets_from_file(filename);

% Let us now extract a tripplet of datasets
alcl_polar_dataset_re3e6clean     = datasets.A01_data;
alcm_polar_dataset_re3e6clean     = datasets.A07_data;
clcd_polar_dataset_re3e6clean     = datasets.B01_data;

% And get some processed polars out of it
[processed_polars_re3e6clean]     = TR824_reader.process_polar_triplet(alcl_polar_dataset_re3e6clean, alcm_polar_dataset_re3e6clean, clcd_polar_dataset_re3e6clean);




% Now l

% % Plot diagnostics
% figure(1)
% plot(al_alcl_polar_sorted    , cl_alcl_polar_sorted    , 'x-'); hold on;
% plot(al_alcl_polar_monotonous, cl_alcl_polar_monotonous, 'o-'); grid on;
% plot(al_filtered             , cl_filtered             , '*-');

figure(2)
subplot(122)
plot(processed_polars_re3e6clean.al_alcl_polar   , processed_polars_re3e6clean.cl_alcl_polar   , 'x-'); hold on;
plot(processed_polars_re3e6clean.al_alclcd_polar , processed_polars_re3e6clean.cl_alclcd_polar , 'o-'); grid on;
xlabel('\alpha (deg)'); ylabel('C_l');
subplot(121)
plot(processed_polars_re3e6clean.cd_clcd_polar   , processed_polars_re3e6clean.cl_clcd_polar   , 'x-'); hold on;
plot(processed_polars_re3e6clean.cd_alclcd_polar , processed_polars_re3e6clean.cl_alclcd_polar , 'o-'); grid on;
xlabel('C_d')         ; ylabel('C_l');
legend('Sorted Data' , 'Reinterpolated Data', 'Location', 'East')


% Create interpolants 







% 
%Re_dataset_01 = 3e6;
%cl_polar_alpha_data_dataset_01

%dataset_01_alpha_cl = ;
%dataset_01_alpha_for_cl = ;
%dataset_01_alpha_for_cl = ;

