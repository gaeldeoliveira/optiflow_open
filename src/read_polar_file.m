function raw_polar = read_polar_file(file_path, varargin)
% RAW_POLAR = READ_POLAR_FILE(FILE_PATH , ...)
% RAW_POLAR = READ_POLAR_FILE(FILE_PATH , TARGET_APPLICATION)
%
% Reads data from
%

if ~isempty(varargin)
    target_application = varargin{1};
else
    target_application = 'xfoil';
end


% Set data importation parameters according to application

DELIMITER = ' ';

% if strcmp(target_application(1:5) , 'rfoil')
%     HEADERLINES = 13;
% else
%     HEADERLINES = 12;
% end
HEADERLINES = 13;

% 13 for Rfoil, 12 for Xfoil

%% Import the file
% Deprecated check condition
% %file_path
% fid = fopen(file_path);
% 
% % disp(file_path);
% if fid>0
%     fclose(fid);
% end

% Try to import data from file, catch error if something goes wrong
% (non-existing or invalid format file)
try
    imported_data = importdata(file_path, DELIMITER, HEADERLINES);
catch
    imported_data = [];                                     % set as empty for robustness.
    disp(['Could not load a polar file :' file_path])
end

% Check that file was valid
if ~isfield(imported_data, 'textdata')
    % if not, return an empty raw polar structure
    raw_polar.alpha = [];
    raw_polar.cl = [];
    raw_polar.cd = [];
else
    %% Now read and separate the file header correctly
    % First get the right header line
    remainder = imported_data.textdata{end-1,1};
    
    % Initialize Loop to separate headers (importdata does not do it right)
    headers = cell(0,0);
    n_word = 1;
    
    while (any(remainder))
        % Chop non-whitespace characters
        [chopped,remainder] = strtok(remainder); %#ok<STTOK>
        % Store Chopped into headers cell array, making everything lower case
        headers{n_word} = lower(chopped);
        n_word = n_word+1;
    end
    
    
    %% Now assign data to raw polar fields
    
    % First take alpha column for reordering in alpha
    [ ~  , sort_index ] = sort(imported_data.data(: , 1));
    
    for n_column = 1:(length(headers)-1)
        % Filter out header names for data consistency between xfoil and rfoil
        if  strcmp(headers{n_column} , 'top_xtr')
            headers{n_column} = 's_xtr';
        end
        if  strcmp(headers{n_column} , 'bot_xtr')
            headers{n_column} = 'p_xtr';
        end
        if  strcmp(headers{n_column} , 're(cl)')
            headers{n_column} = 're_cl';
        end
        
        % Assign data
        raw_polar.(headers{n_column}) = imported_data.data(sort_index , n_column);
    end
    
    % Typically, for Rfoil:
    %   Column 1 : alpha
    %   Column 2 : cl
    %   Column 3 : cd
    %   Column 4 : Re(Cl)
    %   Column 5 : cm
    %   Column 6 : s_xtr
    %   Column 7 : p_xtr
    %   Column 8 : cdp
    % Xfoil only has 7 columns (no Re(Cl)) and order changes, with Cdp on
    % column 4. Also, in Xfoil, s_xtr/p_str corresponds to top_str/bot_str
    
    
    
    
end


%% Old Code Bits
% raw_polar.alpha = imported_data.data(: , 1);                % Extract Angle of attack column
% raw_polar.cl    = imported_data.data(: , 2);                % Extract Cl column
% raw_polar.cd    = imported_data.data(: , 3);                % Extract Cd column
% raw_polar.cm    = imported_data.data(: , 5);                % Extract Cd column
% raw_polar.s_xtr = imported_data.data(: , 6);
% raw_polar.p_xtr = imported_data.data(: , 7);
% raw_polar.cdp   = imported_data.data(: , 8);




% fid = fopen(file_path);  % Open file
%
% if fid > 0
%
%     for k=1:13                              % Skip file header
%     % 13 for Rfoilsuc
%     % 12 for Xfoil
%         fgetl(fid);
%     end
%
%     data =fscanf(fid,'%f');                 % Import data from file into a sequential vector
%     fclose(fid);                            % Close file
%
%     data = transpose(reshape(data, 7, length(data) / 7));   % Reshape data for easy manipulation (7 is due to the number of columns in the data)
%
%     raw_polar.alpha = data(: , 1);                % Extract Angle of attack column
%     raw_polar.cl    = data(: , 2);                % Extract Cl column
%     raw_polar.cd    = data(: , 3);                % Extract Cd column
%     raw_polar.cm    = data(: , 5);                % Extract Cd column
%
% else
%     raw_polar = struct();
% end


end

