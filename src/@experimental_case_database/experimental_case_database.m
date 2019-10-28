classdef experimental_case_database < handle
    %EXPERIMENTAL CASE DATABASE is a class designed to organize a
    % collection of experimental cases, which contain results of wind 
    % tunnel tests on airfoils.
    %
    %   Most fields should be considered as private, but they are left
    %   public for making consultations easier during development phases
    
    properties
        % Fields that should not be written by anyone by the class itself 
        uuid                % [uuid ] - Unique identifier of class         
        EC_cell             % [EC   ] - Cell array of experimental_cases        
        valid_case_index    % [int  ] - Int. array of database indices that refer to valid experimental cases ()
        
        % Helper objects
        SC                  % System context object
        SD                  % Shape definition object
        SF                  % Shape fit object
        
        % Fields that can be taylored externally
        db_folder   = ['./data/experimental/'    ];
        
    end
    
    methods
        function EDB = experimental_case_database()
            % Constructor for the experimental case database, used on creation. 
            % Not yet decided if it will also be used on Construct an instance of this class
            EDB.uuid = EDB.make_a_uuid();
            % Make the helper objects
            EDB.make_helper_objects()
        end
        
        function make_helper_objects(EDB)
            % Create System Context Object 
            EDB.SC = system_context(); 
            % Set for sequential running and set context
            EDB.SC.N_cores = 1; EDB.SC.set_context();
            
            %% Instanciate geometric manipulations classes
            % Parametrization objects for upper and lower side
            p_upper = parametrization( 'cst_upper'  , 12);
            p_lower = parametrization( 'cst_lower'  , 12);
            
            % Fix number of non-CST-shape parameters (1 for TE thickness)
            N_dummy_parameters = 1;
            
            % Shape Definition Objects using previously defined parametrizations
            EDB.SD    = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'cst88');
            
            % Make Shape Fit Object
            EDB.SF = shape_fit_cst(EDB.SD, 'cst88');
        end
                
        function result_flag = add_experimental_case(EDB, EC)
            % Check that new case is not in database yet
            index_list = EDB.find_cases_by_uuid(EC.uuid);
            % Add it if it is not in there
            if isempty(index_list)
                EDB.EC_cell{length(EDB.EC_cell) + 1} = EC;
                result_flag = true;
                lg.msg(EDB, ['Experimental Case added, uuid=' EC.uuid]);
            else
                lg.msg(EDB, ['Experimental Case already exists, uuid=' EC.uuid]);
                result_flag = false;
            end
        end
        
        function result_flag = remove_experimental_case(EC)
            % Find case
            index_list = EDB.find_cases_by_uuid(EC.uuid);
            % Remove if it was found
            if not(isempty(index_list))
                n_index = 1;
                EDB.EC_cell([1:(index_list(n_index)-1) , (index_list(n_index)+1):end]);
                result_flag = true;
                lg.msg(EDB, ['Experimental Case removed, uuid=' EC.uuid]);
            else
                lg.msg(EDB, ['Experimental Case could not be found, uuid=' EC.uuid]);
                result_flag = false;
            end
            
            
        end
        
        function index_list = find_cases_by_uuid(EDB, uuid)
            % Returns a cell array 
            index_list = [];
            for n_EC = 1:length(EDB.EC_cell)
                if isequal(uuid, EDB.EC_cell{n_EC}.uuid)
                    index_list = [index_list, n_EC]; %#ok<AGROW>
                end
            end
        end
        
        function index_list = find_cases_by_airfoil_name(EDB, airfoil_name)
            % Returns a cell array 
            index_list = [];
            for n_EC = 1:length(EDB.EC_cell)
                if isequal(airfoil_name, EDB.EC_cell{n_EC}.airfoil_name)
                    index_list = [index_list, n_EC]; %#ok<AGROW>
                end
            end
        end
        
        function index_list = find_cases_by_IDkind(EDB, IDkind)
            % Returns a cell array 
            index_list = [];
            for n_EC = 1:length(EDB.EC_cell)
                if isequal(IDkind, EDB.EC_cell{n_EC}.IDkind)
                    index_list = [index_list, n_EC]; %#ok<AGROW>
                end
            end
        end
        
        function index_list = find_cases_by_IDkind_on_dbsubset(EDB, IDkind, index_dbsubset)
            % Returns a cell array 
            index_list = [];
            for n_in_dbsubset = 1:length(index_dbsubset)
                n_EC = index_dbsubset(n_in_dbsubset);
                %for n_EC = 1:length(EDB.EC_cell)
                if isequal(IDkind, EDB.EC_cell{n_EC}.IDkind)
                    index_list = [index_list, n_EC]; %#ok<AGROW>
                end
                %end
            end
            
        end
        
        
        function flush_db_snapshot_to_disk(EDB) 
            % Genereate database filename
            db_filename = ['EDB_' , EDB.uuid , '.mat'];
            % Save database to folder
            save([EDB.db_folder , db_filename],'EDB')
        end
    end
    
    methods(Static)
        function uuid = make_a_uuid()
            uuid = char(java.util.UUID.randomUUID);
        end
    end
    
end


