% Script for Reading a Polar File from TR824
% 
%   alfa5 tested manually on 131.ALL and 132.ALL and exhaustively for
%   chapters A and B. Chapter has slightly different filename conventions
%   but we will not use that data explicitly.

% Set Name of File to Read
filename = './data/experimental/NACA/TR824-Digitized/132.ALL';

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

% Close file
fclose(fid);

% Make a consistency check on data set numbers
if (dataset_A01_number ==  1); disp('PASS: dataset_A01 number'); else; disp('FAIL: dataset_A01 number'); end
if (dataset_A02_number ==  2); disp('PASS: dataset_A02 number'); else; disp('FAIL: dataset_A02 number'); end
if (dataset_A03_number ==  3); disp('PASS: dataset_A03 number'); else; disp('FAIL: dataset_A03 number'); end
if (dataset_A04_number ==  4); disp('PASS: dataset_A04 number'); else; disp('FAIL: dataset_A04 number'); end
if (dataset_A05_number ==  5); disp('PASS: dataset_A05 number'); else; disp('FAIL: dataset_A05 number'); end
if (dataset_A06_number ==  6); disp('PASS: dataset_A06 number'); else; disp('FAIL: dataset_A06 number'); end
if (dataset_A07_number ==  7); disp('PASS: dataset_A07 number'); else; disp('FAIL: dataset_A07 number'); end
if (dataset_A08_number ==  8); disp('PASS: dataset_A08 number'); else; disp('FAIL: dataset_A08 number'); end
if (dataset_A09_number ==  9); disp('PASS: dataset_A09 number'); else; disp('FAIL: dataset_A09 number'); end
if (dataset_A10_number == 10); disp('PASS: dataset_A10 number'); else; disp('FAIL: dataset_A10 number'); end
if (dataset_A11_number == 11); disp('PASS: dataset_A11 number'); else; disp('FAIL: dataset_A11 number'); end
if (dataset_A12_number == 12); disp('PASS: dataset_A12 number'); else; disp('FAIL: dataset_A12 number'); end


if (dataset_B01_number ==  1); disp('PASS: dataset_B01 number'); else; disp('FAIL: dataset_B01 number'); end
if (dataset_B02_number ==  2); disp('PASS: dataset_B02 number'); else; disp('FAIL: dataset_B02 number'); end
if (dataset_B03_number ==  3); disp('PASS: dataset_B03 number'); else; disp('FAIL: dataset_B03 number'); end
if (dataset_B04_number ==  4); disp('PASS: dataset_B04 number'); else; disp('FAIL: dataset_B04 number'); end
if (dataset_B05_number ==  5); disp('PASS: dataset_B05 number'); else; disp('FAIL: dataset_B05 number'); end
if (dataset_B06_number ==  6); disp('PASS: dataset_B06 number'); else; disp('FAIL: dataset_B06 number'); end
if (dataset_B07_number ==  7); disp('PASS: dataset_B07 number'); else; disp('FAIL: dataset_B07 number'); end
if (dataset_B08_number ==  8); disp('PASS: dataset_B08 number'); else; disp('FAIL: dataset_B08 number'); end
if (dataset_B09_number ==  9); disp('PASS: dataset_B09 number'); else; disp('FAIL: dataset_B09 number'); end
if (dataset_B10_number == 10); disp('PASS: dataset_B10 number'); else; disp('FAIL: dataset_B10 number'); end
if (dataset_B11_number == 11); disp('PASS: dataset_B11 number'); else; disp('FAIL: dataset_B11 number'); end
if (dataset_B12_number == 12); disp('PASS: dataset_B12 number'); else; disp('FAIL: dataset_B12 number'); end

% 
Re_dataset_01 = 3e6;
cl_polar_alpha_data_dataset_01

dataset_01_alpha_cl = ;
dataset_01_alpha_for_cl = ;
dataset_01_alpha_for_cl = ;

