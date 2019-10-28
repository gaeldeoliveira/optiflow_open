classdef TR824_datasets
    %TR824_DATASETS is a simple class for holding data. For
    %now, we just use it as structure template (non-handle, without methods)
    
    properties
        % Structure datasets has fields:
           % -> A datasets (2 column float arrays)
           A01_data     % alpha, Cl at Re=3e6, clean configuration
           A02_data     % alpha, Cl at Re=6e6, clean configuration
           A03_data     % alpha, Cl at Re=9e6, clean configuration
           A04_data     % alpha, Cl at Re=6e6, rough configuration
           A05_data     % alpha, Cl at Re=6e6, clean configuration, split flap
           A06_data     % alpha, Cl at Re=6e6, rough configuration, split flap
           A07_data     % alpha, Cm at Re=3e6, clean configuration, around c/4 point
           A08_data     % alpha, Cm at Re=6e6, clean configuration, around c/4 point
           A09_data     % alpha, Cm at Re=9e6, clean configuration, around c/4 point
           A10_data     % alpha, Cm at Re=6e6, rough configuration, around c/4 point
           A11_data     % alpha, Cm at Re=6e6, clean configuration, around c/4 point, split flap
           A12_data     % alpha, Cm at Re=6e6, rough configuration, around c/4 point, split flap
           % -> B datasets (2 column float arrays)
           B01_data     % Cl, Cd at Re=3e6, clean configuration
           B02_data     % Cl, Cd at Re=6e6, clean configuration
           B03_data     % Cl, Cd at Re=9e6, clean configuration
           B04_data     % Cl, Cd at Re=6e6, rough configuration
           B05_data     % Cl, Cd at Re=6e6, clean configuration, split flap
           B06_data     % Cl, Cd at Re=6e6, rough configuration, split flap
           B07_data     % Cl, Cm at Re=3e6, clean configuration, around mac point
           B08_data     % Cl, Cm at Re=6e6, clean configuration, around mac point
           B09_data     % Cl, Cm at Re=9e6, clean configuration, around mac point
           B10_data     % Cl, Cm at Re=6e6, rough configuration, around mac point
           B11_data     % Cl, Cm at Re=6e6, clean configuration, around mac point, split flap
           B12_data     % Cl, Cm at Re=6e6, rough configuration, around mac point, split flap
           % -> Information about data source
           filename     % Name of filename from which TR824 data was read
           source='TR824' % Type identifier
    end
    
    methods
        function datasets = TR824_datasets()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end

