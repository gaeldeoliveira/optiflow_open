% Script for Reading a Polar File from TR824
%
%

% Set Name of File to Read
filename = '/home/gael/Desktop/Repositories/optiflow_release/data/experimental/NACA/TR824-Digitized/132.ALL';

% Open the polar file for reading ()
fid = fopen(filename, 'r');

% Read first header line
line_first_header               = fgetl(fid);
% Read first dataset separator
line_first_separator            = fgetl(fid);

% % Start reading dataset
% Read first dataset
% [dataset_number, dataset_data] = TR824_reader.read_dataset_data(fid);

% Read alpha vs cl (1-6) and cm (7-12) datasets
[dataset_A01_number, dataset_A01_data] = TR824_reader.read_dataset_data(fid);
[dataset_A02_number, dataset_A02_data] = TR824_reader.read_dataset_data(fid);
[dataset_A03_number, dataset_A03_data] = TR824_reader.read_dataset_data(fid);
[dataset_A04_number, dataset_A04_data] = TR824_reader.read_dataset_data(fid);
[dataset_A05_number, dataset_A05_data] = TR824_reader.read_dataset_data(fid);
[dataset_A06_number, dataset_A06_data] = TR824_reader.read_dataset_data(fid);
[dataset_A07_number, dataset_A07_data] = TR824_reader.read_dataset_data(fid);
[dataset_A08_number, dataset_A08_data] = TR824_reader.read_dataset_data(fid);
[dataset_A09_number, dataset_A09_data] = TR824_reader.read_dataset_data(fid);
[dataset_A10_number, dataset_A10_data] = TR824_reader.read_dataset_data(fid);
[dataset_A11_number, dataset_A11_data] = TR824_reader.read_dataset_data(fid);
[dataset_A12_number, dataset_A12_data] = TR824_reader.read_dataset_data(fid);

% Read first header line
line_second_header              = fgetl(fid);
% Read first dataset separator
line_second_separator            = fgetl(fid);
% Read Cl, Cd (1-6) and Cl, Cm (7-12) datasets
[dataset_B01_number, dataset_B01_data] = TR824_reader.read_dataset_data(fid);
[dataset_B02_number, dataset_B02_data] = TR824_reader.read_dataset_data(fid);
[dataset_B03_number, dataset_B03_data] = TR824_reader.read_dataset_data(fid);
[dataset_B04_number, dataset_B04_data] = TR824_reader.read_dataset_data(fid);
[dataset_B05_number, dataset_B05_data] = TR824_reader.read_dataset_data(fid);
[dataset_B06_number, dataset_B06_data] = TR824_reader.read_dataset_data(fid);
[dataset_B07_number, dataset_B07_data] = TR824_reader.read_dataset_data(fid);
[dataset_B08_number, dataset_B08_data] = TR824_reader.read_dataset_data(fid);
[dataset_B09_number, dataset_B09_data] = TR824_reader.read_dataset_data(fid);
[dataset_B10_number, dataset_B10_data] = TR824_reader.read_dataset_data(fid);
[dataset_B11_number, dataset_B11_data] = TR824_reader.read_dataset_data(fid);
[dataset_B12_number, dataset_B12_data] = TR824_reader.read_dataset_data(fid);


% Read second header line
line_third_header              = fgetl(fid);
% Read second dataset separator
line_third_separator           = fgetl(fid);
% Read alfa, Cl, Cd (1-6) and alfa, Cl, Cm (7-12) datasets (not sure this
% is reliable! Data interpolation not so good, and file format is
% different, one less line on empty datasets)%
[dataset_C01_number, dataset_C01_data] = TR824_reader.read_dataset_data(fid);
[dataset_C02_number, dataset_C02_data] = TR824_reader.read_dataset_data(fid);
[dataset_C03_number, dataset_C03_data] = TR824_reader.read_dataset_data(fid);
[dataset_C04_number, dataset_C04_data] = TR824_reader.read_dataset_data(fid);
[dataset_C05_number, dataset_C05_data] = TR824_reader.read_dataset_data(fid);
[dataset_C06_number, dataset_C06_data] = TR824_reader.read_dataset_data(fid);
[dataset_C07_number, dataset_C07_data] = TR824_reader.read_dataset_data(fid);
[dataset_C08_number, dataset_C08_data] = TR824_reader.read_dataset_data(fid);
[dataset_C09_number, dataset_C09_data] = TR824_reader.read_dataset_data(fid);
[dataset_C10_number, dataset_C10_data] = TR824_reader.read_dataset_data(fid);
[dataset_C11_number, dataset_C11_data] = TR824_reader.read_dataset_data(fid);
[dataset_C12_number, dataset_C12_data] = TR824_reader.read_dataset_data(fid);


% Close file
fclose(fid);

