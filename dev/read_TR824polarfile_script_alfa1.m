% Script for Reading a Polar File from TR824
%
%

% Set Name of File to Read
filename = '/home/gael/Desktop/Repositories/optiflow_release/data/experimental/NACA/TR824-Digitized/131.ALL';

% Open the polar file for reading ()
fid = fopen(filename, 'r')

% Read first line
line0 = fgetl(fid)
% % Move on to next line
% Read dataset separator
line1 = fgetl(fid)

% % Start reading dataset
% Initialize stuff

% Dataset number
line_dataset_number = fgetl(fid)
dataset_number      = sscanf(line_dataset_number, '%f')

% Data line
line_dataset_data   = fgetl(fid)
dataset_data        = sscanf(line_dataset_data, '%f', 2)


function [line, data] = read_dataset_data_line()
    line = fgetl(fid)
    data = sscanf(line_dataset_data, '%f', 2)
end