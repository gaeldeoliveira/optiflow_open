function ZPGyaw0baseo = read_header_from_dat_file(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   ZPGYAW0BASEO = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   ZPGYAW0BASEO = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   ZPGyaw0baseo = importfile('ZPG_yaw0_base_o.dat', 1, 3);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/01/27 20:18:41

%% Initialize variables.
delimiter = {',','='};
if nargin<=2
    startRow = 1;
    endRow = 3;
end

%% Format string for each line of text:
%   column2: text (%s)
%	column3: text (%s)
%   column4: text (%s)
%	column5: text (%s)
%   column6: text (%s)
%	column7: text (%s)
%   column8: text (%s)
%	column9: text (%s)
%   column10: text (%s)
%	column11: text (%s)
%   column12: text (%s)
%	column13: text (%s)
%   column14: text (%s)
%	column15: text (%s)
%   column16: text (%s)
%	column17: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
ZPGyaw0baseo = [dataArray{1:end-1}];
