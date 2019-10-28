classdef experimental_case_bl_run < handle
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % For templating
        uuid                        % [uuid  ] - UUID of experimental case 
        EDB                         % [class ] - Experimental case database to which this case belongs
        EDB_uuid                    % [uuid  ] - UUID of experimental case database
        
        % Specific
        data_source                 % [string] - String with source of data (for UI search and display purposes only)
    end
    
    properties (SetAccess= private)
        % For templating
        IDkind = 'bl_run'           % [string] - Protected string describes the type of case this describes
        
        % Specific: Structure data stores
        bl_run                      % [obj.  ] - Description (templated) of experimental bl run
    end
    
    methods
        function EC = experimental_case_bl_run(EDB)
            % EXPERIMENTAL_CASE_BL_RUN constructor method. (Part of this
            % procedure should be templated at a later stage! To ensure
            % consistency accross experimental_cases of all kinds!)
            %   Instances of the experimental case class are constructed
            %   without providing data. Instead, they receive a hook to an
            %   an experimental case database object, which guides them on
            %   how to process data in a standardized way. 
            %
            %   Even though the EDB is not a parent in the formal sense of
            %   the world, it can be thought of as such because it governs
            %   the behavior of the EC to a large extent.
            %   
            %   Data is stored in structured templates.
            
            % Store hook to EDB object 
            EC.EDB = EDB;
            EC.EDB_uuid = EC.EDB.uuid;
            % Make uuid of experimental case
            EC.uuid = EC.EDB.make_a_uuid();
            % Initialize appropriate flags to boolean 
        end
        
        function load_bl_run(EC, bl_run)
            % Store Object
            EC.bl_run = bl_run;
            % Add data source display
            EC.data_source = bl_run.source;
        end
        
        function standardize_bl_run(EC) %#ok<MANU>
            % Write this function in due time!
        end
    end
end

