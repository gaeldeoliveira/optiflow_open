classdef shape_fit < handle
    %SHAPE_DEFINITION This class stores, aquires and converts airfoil shape
    % definitions. It is a superclass for parametrization specific fitting
    % classes
     
    properties 
        parameters
        shape_definition_handle        
        reference_source        
        imported_coordinates        
        name
        id
%       parameters_down
%       parameters_suction

    end
    
    methods      
        function import_coordinates(SD, fileToRead)
            % Import coordinates from airfoil file, treated as Lednicer
            % (standard Xfoil/Rfoil) for all file extensions except '.gnu'
            % which indicates a Selid file as outputed from NACA456
            % shape generator (PDAS version, Carmichael)
            
            if strcmp(fileToRead(end-3:end), '.gnu')
                SD.import_selig_coordinates(fileToRead)
            else
                SD.import_lednicer_coordinates(fileToRead)
            end
        end
        
        function import_lednicer_coordinates(SD, fileToRead)
            % Import coordinates from headerfree Lednicer airfoil file 
            
            % fileToRead1 = '/home/gael.deoliveira/Documents/MATLAB/actiairfoil_0.1alpha/data/airfoils/du97W300.txt';
            % disp(fileToRead)
            newData1 = importdata(fileToRead);
            
            % Convert to structure for cases when file had no header.
            if ~isstruct(newData1)
                data = newData1;
                newData1 = struct;
                newData1.data = data;
            end
            % Store imported data as it is in file (for verification
            % purposes)
            SD.imported_coordinates.raw.tx_coordinates = newData1.data(:,1);
            SD.imported_coordinates.raw.tz_coordinates = newData1.data(:,2);
            
            % Preprocess Data to ensure consistent position and chord
            SD.preprocess_coordinates();            
        end
        
        
        function import_selig_coordinates(SD, fileToRead)
            % Import coordinates from Selig airfoil file, as outputed from
            % NACA456 shape generator. LE assumed to be exactly (in binary,
            % not EPS sense!) at (0,0) and supplied on both sides. Headers
            % are currently unsupported.
            %
            % This file is working a bit as a workaround, in the sense that
            % 
                         
            % fileToRead1 = '/home/gael.deoliveira/Documents/MATLAB/actiairfoil_0.1alpha/data/airfoils/du97W300.txt';
            % disp(fileToRead)
            newData1 = importdata(fileToRead);
            
            % Extract x coordinates (for both sides, from LE to TE for both sides)
            tx_selig = newData1(:,1);
            tz_selig = newData1(:,2);
            
            % Find indices of leading edge for each side (points where both
            % x and y are zero)
            le_indices   = find(and(tx_selig  == 0, tz_selig == 0));
            le_top_index = le_indices(1);
            le_bot_index = le_indices(2);
            
            % Make upper selig coordinates
            tx_top = tx_selig(le_top_index:(le_bot_index-1)); 
            tz_top = tz_selig(le_top_index:(le_bot_index-1));
            % Make lower selig coordinates
            tx_bot = tx_selig(le_bot_index:end);
            tz_bot = tz_selig(le_bot_index:end);

            % Collate into a lednicer set
            raw_lednicer_tx_coordinates = [flipud(tx_top) ; tx_bot(2:end)];
            raw_lednicer_tz_coordinates = [flipud(tz_top) ; tz_bot(2:end)];

            % Store imported data as it is in file (for verification
            % purposes)
            SD.imported_coordinates.raw.tx_coordinates = raw_lednicer_tx_coordinates;
            SD.imported_coordinates.raw.tz_coordinates = raw_lednicer_tz_coordinates;
            
            % Preprocess Data to ensure consistent position and chord
            SD.preprocess_coordinates(); 
        end
        
        
        function import_ISES_coordinates(SD, fileToRead)
            % fileToRead1 = '/home/gael.deoliveira/Documents/MATLAB/actiairfoil_0.1alpha/data/airfoils/du97W300.txt';
            newData1 = importdata(fileToRead);
            
            % Store imported data as it is in file (for verification
            % purposes)
            SD.imported_coordinates.raw.tx_coordinates = newData1.data(3:end,1);
            SD.imported_coordinates.raw.tz_coordinates = newData1.data(3:end,2);
            
            % Preprocess Data to ensure consistent position and chord
            SD.preprocess_coordinates();
        end
        function preprocess_coordinates(SD)
%             % Preprocess data to ensure that leading edge is at tx=0 tz=0 and
%             % trailing edge at tx = 1
%             %    First determine translation and scaling constants
%             translation_factor = min(SD.imported_coordinates.raw.tx_coordinates);
%             scale_factor = max(SD.imported_coordinates.raw.tx_coordinates) - min(SD.imported_coordinates.raw.tx_coordinates);
%             %    Now apply translation on x coordinate, followed by scaling
%             %    on both coordinates
%             SD.imported_coordinates.adimensionalized.tx_coordinates = (SD.imported_coordinates.raw.tx_coordinates - translation_factor) / scale_factor;
%             SD.imported_coordinates.adimensionalized.tz_coordinates = SD.imported_coordinates.raw.tz_coordinates / scale_factor;
%             
%             % Now add linear correction factor to ensure that z_te is the
%             % same on both sides

            % Find TE coordinates
            TE_x = max(SD.imported_coordinates.raw.tx_coordinates(1) , SD.imported_coordinates.raw.tx_coordinates(end));    % Handle case in which upper and lower TE do not have the same coordinate
            TE_z = 0.5 * (SD.imported_coordinates.raw.tz_coordinates(1) + SD.imported_coordinates.raw.tz_coordinates(end));  % Define TE such that upper and lower TE bluntness are the same

            % Now translate TE so that TE is at (tx,tz) = (1,0)
            tx = SD.imported_coordinates.raw.tx_coordinates - TE_x + 1;
            tz = SD.imported_coordinates.raw.tz_coordinates - TE_z;

            % Now identify LE
            % Find Point that is most distant from T.E.
            [original_chord, mindex] = max(( tx - 1).^2 + (tz - 0).^2);

            % Find angle by which airfoil must be rotated to get plane chord
            % (approximate method valid only for small angles)
            th = 2 * asin(tz(mindex) / (2* original_chord));

            % Now rotate coordinates of airfoil
            [tx, tz, ~] = rotatearbitrary(  tx, tz, 0 , 1 , 0 , 0 , 0,0,1, th);

            % Scale airfoil and translate airwfoil coordinates to final position
            scale_factor = max(tx) - min(tx);
            translation_factor = min(tx);
        
            tx = (tx-translation_factor) / scale_factor; % Preserve TE position
            tz = tz / scale_factor;

            % Now store!
            SD.imported_coordinates.adimensionalized.tx_coordinates = tx;
            SD.imported_coordinates.adimensionalized.tz_coordinates = tz;
            
            % Add storage of modifications from raw coordinates
            SD.imported_coordinates.processing_parameters.theta              = th;
            SD.imported_coordinates.processing_parameters.scale_factor       = scale_factor;
            SD.imported_coordinates.processing_parameters.translation_factor = translation_factor;            
        end

    end
    % Restons Amants - Maxime le Forestier
end

