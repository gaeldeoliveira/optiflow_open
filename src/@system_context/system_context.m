classdef system_context < handle
    % SYSTEM_CONTEXT class stores data about current system, in order to 
    %   facilitate multiOS operation (even though windows support is not
    %   very solid! and linux is always a better choice). Viva o pinguim!
    
    properties
        base_dir = [];
        src_subdir = ['src' filesep()];     % Generate OS adap
        airfoil_subdir = ['data' filesep() 'airfoils' filesep()];
        feasible_airfoil_subdir = ['data' filesep() 'feasible_airfoils' filesep()];
        tmp_subdir = [ 'tmp' filesep() datestr(now, 30) filesep()];      % Set time on construction to avoid file mixing between cases!
        apps_subdir = ['apps' filesep()];
        core_subdirs = [];
        core_dir_prefix = 'core';
        user_src_subdir = ['user_src' filesep()];
        results_subdir = ['results' filesep()];
        
        airfoil_filename = 'xfin';
        polar_filename = 'plr';
        polar_dump_filename = 'dump';
        svbl_filename = 'svbl';
        svcp_filename = 'cpx';
        
        app_extension = '';            % Application Extension, to keep multi-OS app binaries at hand!
        
        hObject = [];                  % Handle to graphical user interface
        
        creation_time_str = datestr(now, 30);
        suc_filename = 'suc';
        
        fs = filesep();                % Folder path separator on current system (default for Linux), changed automatically for Windows
        fs_sprinf_ap                   % Folder path separator appendix for Sprintf MultiOS compatilibity (='' on unix, ='\' on windows) 
        OS = computer();               % Are we on a *NIX or Win Environment? Important for file handling
        N_cores = 6;                   % Number of Cores to be used during calculation
        
    end
    
    methods
        function  SC = system_context() 
            % Find my own path/folder
            p_full = mfilename('fullpath');
            p_file = mfilename();
            SC.base_dir = p_full(1:(end-(length(p_file)*2 + 2 + length(SC.src_subdir))));
            % Set project path as MATLAB's current directory
            cd(SC.base_dir)
            % Set fs_sprinf_ap in function of system
            if ispc
                SC.fs_sprinf_ap = filesep();
            else
                SC.fs_sprinf_ap = '';
            end
        end
        
        function SC = set_context(SC)
            % Add extra directories to path
            
            % Find my own path/folder
            cd(SC.base_dir)
            
            % Make tmp subdirectory
            mkdir(SC.base_dir , SC.tmp_subdir)
            
            SC.core_subdirs = cell(SC.N_cores,1);
            % Make core subdirectories
            for n_core  = 1:SC.N_cores
                mkdir([SC.base_dir , SC.tmp_subdir], [SC.core_dir_prefix  , num2str(n_core)])
                SC.core_subdirs{n_core} = [SC.core_dir_prefix  , num2str(n_core) , SC.fs];
            end
            
            % Add user cost functions to matlab search path
            addpath([ SC.base_dir , SC.user_src_subdir])
            
            % Make extension of rfoil application
            SC.make_app_extension();

            % Set Environment Variables for MAC (taylor /usr/lib to your gfortran
            % installation)
            if ismac
                setenv('DYLD_LIBRARY_PATH' , ['/usr/lib' getenv('DYLD_LIBRARY_PATH')]);
            end
            
            % Disable some usual interpolation warnings to ensure readable
            % output
            warning('off', 'MATLAB:interp1:ppGriddedInterpolant');
            warning('off', 'MATLAB:interp1:UsePCHIP');
        end
        
        function SC = migrate_context(SC, identifier)
            % Find my own path/folder
            p_full = mfilename('fullpath');
            p_file = mfilename();
            SC.base_dir = p_full(1:(end-(length(p_file)*2 + 2 + length(SC.src_subdir))));
            % Set project path as MATLAB's current directory
            cd(SC.base_dir)
            % Edit name of tmp_subdir to something (not yet used and
            % customizable)
            SC.tmp_subdir = [ 'tmp' filesep() datestr(now, 30) identifier filesep()]; 
            % Set (new) context
            SC = SC.set_context();
            % Inform user:
            disp('SC: Migrated SC (and hence case) to:')
            disp([SC.base_dir SC.tmp_subdir])
           
        end
        
        function make_app_extension(SC)
            if ispc 
                SC.app_extension = '.exe';  % Windows App
            end
            if ismac
                SC.app_extension = '.mac';  % Mac OsX App
            end
            if and(isunix, not(ismac))
                SC.app_extension = '.lin';  % Linux App (is unix, but not mac)
            end
        end
    end
    
end

