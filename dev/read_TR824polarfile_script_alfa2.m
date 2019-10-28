% Script for Reading a Polar File from TR824
%
%

% Set Name of File to Read
filename = '/home/gael/Desktop/Repositories/optiflow_release/data/experimental/NACA/TR824-Digitized/131.ALL';

% Open the polar file for reading ()
fid = fopen(filename, 'r');

% Read first line
line0 = fgetl(fid);
% % Move on to next line
% Read dataset separator
line1 = fgetl(fid);

% % Start reading dataset
% Read first dataset
[dataset_number, dataset_data] = TR824_reader.read_dataset_data(fid);

% And do this manually again a few times!

% Close file
fclose(fid);

