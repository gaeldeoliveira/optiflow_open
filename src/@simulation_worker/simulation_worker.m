classdef simulation_worker < handle
    % SIMULATION_WORKER class is a proletarian dedicated to executing all
    % sorts of simulations described in simul_protocol, on airfoils
    % described by a shape definition SD, within a system context SC.
    % Optionally, airfoils from files can be analyzed, in which case SD can
    % be left empty (replaced by [ ]).    
    %
    % Simulation worker has lived through the industrial revolution,
    % whereby it is capable of parallel execution, when the parallelized
    % field is set to a number greater than 0.
    %
    % Before its methods are used the class has to be constructed using its
    % primary method declared as:
    %    function  SW = simulation_worker(name, simul_protocol , SPD , SD, SC, SDD (optional) )
    %
    % Despite being a strong mass, proletarians are neither homogeneous nor
    % anonymous, whereby they  require varied strings for the name field to
    % avoid intereferences.
    % The SDD object link is a varargin and helps the SW modify the shape
    % it is analyzing. See SDD class description for more.
   
    properties
        % Defined on Creation
        name                             % Name of worker (gives name to temporary script)
        simul_protocol                   % Handle to simulation protocol
        SPD                              % Handle to simulation protocol dynamizer
        SD                               % Handle to shape definition 
        SC                               % Handle to system context        
        SDD
        
        multisim_id_value                % Passed to the results structure to identify/order results in multsim cases (with more than one simulation worker!)
        
        % Post creation user options
        app_name = 'rfoilsucplasma';     % Name of application to use for execution
        parallelized = 1;                % Enable parallelized execution of simulations (cost function calls must be vectorized, at GM (@gamultiobj_manager) level, for this to have any effect)
        doubleIO = 1;                    % Preserve double precision in Matlab file IO (only effective when app IOs in double precision as well)
        max_wait = 100;                  % Maximum number of seconds to wait for a block of simulations to be executed
        
        % Internals        
        fid                              % FID for writting BASH script file        
        fig_airfoil  = []                % handle to airfoil figure
        fig_polar    = []                % handle to polar figure
        fig_pressure = []                % handle to airfoil pressure distribution
        fig_active   = true;             % flag to activate/deactivate plotting (true by default for GUI) 
        
    end
    
    methods
        function  SW = simulation_worker(name, simul_protocol , SPD , SD, SC, varargin)
            % Constructor method (adapts to current OS, dirty approach)
            SW.name = name;
            SW.simul_protocol = simul_protocol;
            SW.SPD = SPD;
            SW.SD = SD;
            SW.SC = SC;
            % Add simulation/shape dynamizer if it is present
            if not(isempty(varargin))
                SW.SDD = varargin{1};
            end
        end
        
        function  write_shell_script(SW, varargin)
            %     % A typical script we want:
            %     #!/bin/bash
            %     cd ./tmp
            %     cd ./core1
            %     ../../apps/rbink90DPv2options < ../clean_no_suction &
            %     cd ../core2
            %     ../../apps/rbink90DPv2options < ../clean_no_suction &
            %     wait
            
            % If nothing is said execute a run for each core
            if isempty(varargin)
                N_cores_run = SW.SC.N_cores;                
            else
                N_cores_run = varargin{1};
            end
                
            
            % Dummy constructor method
            % LINUX and MAC OSX version
            if not(ispc)
                % Open File
                SW.fid = fopen([SW.SC.tmp_subdir , SW.name], 'wt');
                
                % Write script header
                % #!/bin/bash
                fprintf(SW.fid, '#!/bin/bash \n');
                fprintf(SW.fid, '# Paralel Airfoil Optimization System - Gael de Oliveira \n');
                fprintf(SW.fid, '# Automatically Generated BASH script by the Simulation Worker Class \n \n');
                
                % Enter tmp subdir
                % fprintf(SW.fid, [ 'cd ./'  , SW.SC.tmp_subdir , '  \n']);
                
                % Write each core's block
                for n_core = 1:N_cores_run
                    SW.write_shell_script_call_block(n_core)
                end
                
                % Now wait for termination of all processes from our shell (se
                % we know when we can read the polars in matlab
                fprintf(SW.fid, 'wait   \n');
                
                % We're done with the command file! Close it!
                fclose(SW.fid);
                %             if SW.fid  > 0
                %                 % Only if we are writting to a succesfuly opened file
                %                 fclose(SW.fid);
                %             end
            end
            
            % WINDDOWS VERSION
            if ispc
                for n_core = 1:N_cores_run     
                    % Open File
                    SW.fid = fopen([SW.SC.tmp_subdir , SW.name '_' num2str(n_core) '.bat'], 'wt');
                    % Write Bal Script
                    SW.write_shell_script_call_block(n_core)
                    % We're done with the command file! Close it!
                    fclose(SW.fid);
                end   
            end

        end
        
        function  write_shell_script_call_block(SW , n_core)           
            
            %%% LINUX / MAC OSX Version            
            if not(ispc)
                % Enter Right Core Directory % Deprecated for MultiOS
                % compatibility, because windows shell does not have the 
                % point . and .. operators.
                % fprintf(SW.fid, [ 'cd ./'  , SW.SC.core_subdirs{n_core} , '  \n']);                           % Old Version for single OS (Linux) only
                
                % Prepare Application call line
                % app_str   = [   './' , SW.SC.apps_subdir , SW.app_name];                                      % Old Version for single OS (Linux) only
                app_str = [ './' SW.app_name SW.SC.app_extension];                                              % New MultiOS version

                % Write Tunnel String
                % tunel_str = [   '' , SW.simul_protocol.name];                                                 % Older version for static simulation parameters
                % tunel_str = SW.simul_protocol.name;                                                           % Old version for dynamic simulation parameters, but single OS (Linux)
                tunel_str = ['./' SW.SC.tmp_subdir SW.SC.core_subdirs{n_core} SW.simul_protocol.name];  % New MultiOS version

                % Write dummy_out file string to channel stdout
                dummy_str = ['./' SW.SC.tmp_subdir SW.SC.core_subdirs{n_core}  'dummy_out'];                         % Added for New MultiOS version
                if (SW.parallelized > 0)
                    % app_call_str = [ app_str   , '  <  ' , tunel_str ,  '  >  dummy_out' , ' & \n'];          % Old SingleOS (Linux) version
                    app_call_str = [ app_str   , '  <  ' , tunel_str ,  '  >  ' dummy_str , ' & \n'];          % New MultiOS version
                    % It is necessary to have a file to tunnel stdout file to
                    % avoid matlab mixing up the file IO of the two processes
                    % and other problems i did not understand fully!
                else
                    app_call_str = [ app_str   , '  <  ' , tunel_str ,  ' \n'];
                end
                % Write Application Call line
                fprintf(SW.fid, app_call_str );
                
                % Leave Core Directory
                % fprintf(SW.fid, 'cd ../   \n \n');    % Deprecated in new MultiOS version  
                
            end
            
            %%% WINDOWS Version
            if ispc
                % Prepare Application call line
                app_str   = [ SW.app_name  SW.SC.app_extension];
                
                % It is necessary to have a file to tunnel stdout file to
                % avoid matlab mixing up the file IO of the two processes
                % and other problems i did not understand fully!
                tunel_str = SW.simul_protocol.name;              % New version for dynamic simulation parameters
                if (SW.parallelized > 0)
                    % Duplicate backslashes to avoid ambiguity with of \ with sprintf() control characters (e.g. \n) 
                    app_call_str = [ app_str   , '<'  SW.SC.tmp_subdir  '\' SW.SC.core_subdirs{n_core} , '\' , tunel_str ,  '  > ' SW.SC.tmp_subdir '\' SW.SC.core_subdirs{n_core} '\dummy_out \n'];
                else
                    % Duplicate backslashes to avoid ambiguity with of \ with sprintf() control characters (e.g. \n)
                    app_call_str = [ app_str   , '<'  SW.SC.tmp_subdir  '\' tunel_str ,  '  >'  SW.SC.tmp_subdir '\dummy_out \n '];
                end
                
                % Also prepare semaphor file pathline (the semaphor is used
                % by SW.execute_shell_script() to check for the completion
                % of the batch process (ah... windows....)
                
                
                % Write Application Call line
                fprintf(SW.fid, app_call_str );
                % Write semaphor file creation call (heterodox use of dos shell copy command)
                % Duplicate backslashes to avoid ambiguity with of \ with sprintf() control characters (e.g. \n)
                fprintf(SW.fid, ['copy nul >'  SW.SC.tmp_subdir '\' SW.SC.core_subdirs{n_core}  '\myfile1 \n'] );
                % Write batch file exit clause
                fprintf(SW.fid, 'exit  ' );
            end
            
        end
        
        function  execute_shell_script(SW,varargin)
            
            %%% LINUX / MAC OS X Version
            if not(ispc)
                % Force closing of shell file (correct matlab bug of not releasing file handle)
                status = 1;
                while status ~= 0
                    status = fclose('all');
                end
                
                % Set our script permissions so it is executable
                system(['chmod 777  ' , SW.SC.tmp_subdir, SW.name]);
                % Set environment variables to avoid confusion for Fortran
                % IO routines, Matlab disrupts the Fortran LUs (Logical
                % Unit) numbering, but this did not use to bea problem with 
                % relaxed standard-abidance compilers. It becomes one with 
                % recent gfortran implementations that seem commited to
                % crash all old code's IO.
                setenv('GFORTRAN_STDIN_UNIT'    , '5');         % Equivalent to type 'export GFORTRAN_STDIN_UNIT=5' in BASH
                setenv('GFORTRAN_STDOUT_UNIT'   , '6');
                setenv('GFORTRAN_STDERR_UNIT'   , '0');
                % Call Script for execution
                system(['./',SW.SC.tmp_subdir, SW.name]);
                % Restore environment variables to their usual form for
                % Matlab (Note: shell environment variables propagate from 
                % parent to child but do not have a systemwide impact)
                setenv('GFORTRAN_STDIN_UNIT'    , '-1');        % Equivalent to type 'export GFORTRAN_STDIN_UNIT=-1' in BASH
                setenv('GFORTRAN_STDOUT_UNIT'   , '-1');        % Setting these the standard IO channels to -1 is probably like making a sink, and it therefore makes piping impossible! 
                setenv('GFORTRAN_STDERR_UNIT'   , '-1');
            end
            
            %%% WINDDOWS Version
            if ispc
                % Force closing of shell file (correct matlab bug of not releasing file handle)
                status = 1;
                while status ~= 0
                    status = fclose('all');
                end
                
                % If nothing is said in varargin execute a run for each core
                if isempty(varargin)
                    N_cores_run = SW.SC.N_cores;
                else
                    N_cores_run = varargin{1};
                end
                
                % Dummy constructor method
                for n_core = 1:(N_cores_run-1)
                    % Make hidden file!
                    Temp_batch_file=[SW.SC.tmp_subdir(1:end-1) SW.SC.fs SW.name '_' num2str(n_core)  '.bat &'];
                    % Temp_batch_file=[SW.SC.tmp_subdir(1:end-1) SW.SC.fs SW.name '_' num2str(n_core)  '.bat'];
                    system(Temp_batch_file);
                end
                
%                 parfor n_core = 1:(N_cores_run-1)
%                     % Make hidden file!
%                     % Temp_batch_file=[SW.SC.tmp_subdir(1:end-1) SW.SC.fs SW.name '_' num2str(n_core)  '.bat &'];
%                     Temp_batch_file=[SW.SC.tmp_subdir(1:end-1) SW.SC.fs SW.name '_' num2str(n_core)  '.bat'];
%                     system(Temp_batch_file);
%                 end
                
                
                system([SW.SC.tmp_subdir(1:end-1) SW.SC.fs SW.name '_' num2str(N_cores_run) '.bat']);
                
                % Check for the creation of a file called 'myfile1' whose
                % creation is used as a semaphor to indicate that the
                % batchfiles we called have completed their work (this is
                % part of the compatbility measures which could be made far
                % more solid!)
                for n_core = 1:(N_cores_run)
                    fid_test=-1;
                    while fid_test==-1
                        FILE_POLAR_TEST=[ SW.SC.tmp_subdir SW.SC.core_subdirs{n_core} , 'myfile1'];
                        fid_test_a=fopen(FILE_POLAR_TEST);
                        fid_test_b=fid_test_a;
                        if fid_test_a~=-1
                            fclose(fid_test_a);
                        end
                        pause(.2)
                        fid_test=fid_test_b;
                        pause(.2)
                        disp('Waiting for polar to be calculated... if this lasts too long, stop simulation')
                    end
                end
            end
        end
                
        function  clean_folders(SW, varargin)
            % If nothing is said clean folders for each core
            if isempty(varargin)
                N_cores_run = SW.SC.N_cores;
            else
                N_cores_run = varargin{1};
            end
            
            for n_core = 1:N_cores_run
                % Clean polar file
                delete([SW.SC.tmp_subdir, SW.SC.core_subdirs{n_core} , SW.SC.polar_filename])
                % Clean dump file
                delete([SW.SC.tmp_subdir, SW.SC.core_subdirs{n_core} , SW.SC.polar_dump_filename])
                % Clean dummy_out file
                delete([SW.SC.tmp_subdir, SW.SC.core_subdirs{n_core} , 'dummy_out'])
                % Clean svbl files
                delete([SW.SC.tmp_subdir, SW.SC.core_subdirs{n_core} , SW.SC.svbl_filename , '*'])                
                % WINDOWS Version - Delete extra file for
                if ispc
                    delete([SW.SC.tmp_subdir, SW.SC.core_subdirs{n_core} , 'myfile1'])
                end
            end
        end
        
        function  ap_block_cell = load_polars(SW, varargin)
            % If nothing is said then assume we work on all cores
            
            if isempty(varargin)
                N_cores_run = SW.SC.N_cores;
            else
                N_cores_run = varargin{1};
            end
            
            % Preallocate cell array for aerodynamic polar objects
            ap_block_cell = cell(1, N_cores_run);
            
            for n_core = 1:N_cores_run
                % Load results from polar file
                raw_data = read_polar_file( ...
                    [SW.SC.tmp_subdir, SW.SC.core_subdirs{n_core} , SW.SC.polar_filename],...
                    SW.simul_protocol.target_application);
                
                % Make aerodynamic polar object and store it . Do so in try
                % catch statement to handle file IO errors if R/Xfoil
                % throws strange stuff out
                try
                    ap_block_cell{n_core} = aerodynamic_polar(raw_data);
                catch
                    clear raw_data
                    ap_block_cell{n_core} = NaN;
                    disp('Failed Polar - Could not load raw_data into aerodynamic polar')
                end
                
                % If applicable, load SVBL pressure and boundary layer data 
                
                if strcmp(SW.simul_protocol.save_bl_parameters, 'true')
                    try
                        bldata = read_svbl_files([SW.SC.tmp_subdir, SW.SC.core_subdirs{n_core} , SW.SC.svbl_filename]);
                        ap_block_cell{n_core} = ap_block_cell{n_core}.append_bldata(bldata);
                    catch
                        disp('Failed Polar - Could not read svbl data for appending bldata')
                    end
                end
                
                
            end
        end
        
        function  coordinates_block_cell = write_airfoils_block(SW, x_block_cell)
            % Writes the coordinate files of  a block of airfoils specified 
            % by a cell array of parameters in the form defined by the shape 
            % definition object SD supplied to this worker object
            
            coordinates_block_cell = cell(size(x_block_cell));
            for n_core = 1:length(x_block_cell)
                [tx, tz] = SW.SD.generate_coordinates(...
                    SW.simul_protocol.N_airfoil_points , x_block_cell{n_core});
                
                % Modify airfoil coordinates for analysis according to 
                % simulation protocol specifications
                [analysis_coordinates_cell] = ...
                    SW.simul_protocol.geometry_correction_routine({tx, tz});
                
                tx_analysis = analysis_coordinates_cell{1};                
                tz_analysis = analysis_coordinates_cell{2};
                
                % Airfoil Coordinate definition module
                Xt = [tx_analysis , tz_analysis]; %#ok<NASGU>
                                                
                coordinates_block_cell{n_core}.tx = tx;
                coordinates_block_cell{n_core}.tz = tz;
                coordinates_block_cell{n_core}.tx_analysis = tx_analysis;
                coordinates_block_cell{n_core}.tz_analysis = tz_analysis;
                
                % Build airfoil string file
                airfoil_file = [SW.SC.tmp_subdir, ...
                    SW.SC.core_subdirs{n_core} , SW.SC.airfoil_filename];
                
                % Save to file in right format
                if (SW.doubleIO > 0)
                    save( airfoil_file , 'Xt','-ascii', '-double');
                else
                    save( airfoil_file , 'Xt','-ascii');
                end                               
           end
        end
        
        function  write_simulation_protocols_block(SW, x_block_cell)
            % Writes the (dynamized) simulation protocols into each
            % airfoils folder
            
            % For each airfoil in current execution execution block
            for n_core = 1:length(x_block_cell)                
                % If simulation protocol dynamizer was defined, use it! 
                if ~isempty(SW.SPD)
                    % First dynamize simulation protocol with current airfoil
                    % vector
                    SW.SPD.modify_simulation_protocol(x_block_cell{n_core});                                        
                end
                
                % Then write simulation protocol file using varargin of
                % write_command_file method (by default it writes to the
                % tmp root
                SW.simul_protocol.write_command_file(SW.SC.core_subdirs{n_core});
            end
            
            % We are done!                                                                               
        end            
        
        function  write_simulation_protocols_block_with_polar_conditions(SW, x_block_cell, polar_conditions_block_cell)
            % Writes the (dynamized) simulation protocols into each
            % airfoils folder
            
            % For each airfoil in current execution execution block
            for n_core = 1:length(x_block_cell)                
                % If simulation protocol dynamizer was defined, use it! 
                if ~isempty(SW.SPD)
                    % First dynamize simulation protocol with current airfoil
                    % vector
                    SW.SPD.modify_simulation_protocol(x_block_cell{n_core});                                        
                end
                
                % Edit simulation protocol to mach prescribed polar
                % conditions (this should become a method of the static
                % class experimental_polar_conditions!) 
                SW.simul_protocol.M       = polar_conditions_block_cell{n_core}.M;
                SW.simul_protocol.Re      = polar_conditions_block_cell{n_core}.Re;
                SW.simul_protocol.N_crit  = polar_conditions_block_cell{n_core}.N_crit;
                SW.simul_protocol.xtr_top = polar_conditions_block_cell{n_core}.xtr_top;
                SW.simul_protocol.xtr_bot = polar_conditions_block_cell{n_core}.xtr_bot;
                if (polar_conditions_block_cell{n_core}.INGE == true) % For reverse compatibility!
                    SW.simul_protocol.INGE    = polar_conditions_block_cell{n_core}.INGE;
                end
                
                % Then write simulation protocol file using varargin of
                % write_command_file method (by default it writes to the
                % tmp root
                SW.simul_protocol.write_command_file(SW.SC.core_subdirs{n_core});
            end
            
            % We are done!                                                                               
        end            
        
        function  output_figures(SW , results)
            
%             ap =  results.ap;
%             if ~isempty(SW.fig_polar)
%                ap.plot(SW.fig_polar , 'r');
%             else
%                 SW.fig_polar = ap.plot();
%             end
%             
%             % Output function for monitorization
%             if ~isempty(SW.fig_airfoil)
%                 SW.SD.plot_airfoil(results.x , SW.fig_airfoil);
%             else
%                 SW.fig_airfoil = SW.SD.plot_airfoil(results.x);
%             end
%             
%             % Now single fresh airfoil
%             figure(SW.fig_fresh_airfoil)
%             hold off
%             plot(0,0)
%             SW.SD.plot_airfoil(results.x , SW.fig_fresh_airfoil);

%             % Just call the GUI handler
%             display_results(SW.SC.hObject , results)
            
            % Just call the GUI handler if available
            if ishandle(SW.SC.hObject)
                display_results(SW.SC.hObject , results)
            else
                if (SW.fig_active == true)
                    try
                        if isa(results.ap, 'aerodynamic_polar')
                            if isempty(SW.fig_polar)
                                h4 = figure(4); %#ok<NASGU>
                                SW.fig_polar   = results.ap.plot;
                                legend off;
                                SW.fig_polar(1).Children(1).Color = [0 0.4470 0.7410];
                                SW.fig_polar(2).Children(2).Color = [0 0.4470 0.7410];
                            else
                                SW.fig_polar(1).Children(1).XData = results.ap.cd_alpha([]).*10^4;
                                SW.fig_polar(1).Children(1).YData = results.ap.cl_alpha([]);
                                
                                SW.fig_polar(2).Children(1).XData = results.ap.alpha_cl_max_local_limited;
                                SW.fig_polar(2).Children(1).YData = results.ap.cl_max_local_limited;
                                
                                SW.fig_polar(2).Children(2).XData = results.ap.raw_data.alpha(1):(1/5):results.ap.raw_data.alpha(end);
                                SW.fig_polar(2).Children(2).YData = results.ap.cl_alpha([]);
                            end
                            % %
                            alpha_range = results.ap.alpha_range;
                            ld_alpha = results.ap.cl_alpha(alpha_range) ./ results.ap.cd_alpha(alpha_range);    
                            [ld_max, i_ld_max] = max(ld_alpha); alpha_ld_max = alpha_range(i_ld_max);
                            % %
                            SW.fig_polar(1).Title.String = ['max L/D = ' , num2str(ld_max) , '   at Cl = ' , num2str(alpha_ld_max)];
                            SW.fig_polar(2).Title.String = ['max C_l = ' , num2str(results.ap.cl_max_local_limited)];
                            
                            % % Original Legend
                            SW.fig_polar(1).Title.String = ['max L/D = ' num2str(max(results.ap.cl_alpha([]) ./ results.ap.cd_alpha([])))];
                            SW.fig_polar(2).Title.String = ['max C_l = ' num2str(results.ap.cl_max_local_limited)];
                            
                            if isempty(SW.fig_airfoil)
                                h5 = figure(5); %#ok<NASGU>
                                SW.fig_airfoil = plot(results.coordinates.tx, results.coordinates.tz);
                                grid on; axis equal;
                                xlabel('x/c');ylabel('y/c');
                            else
                                SW.fig_airfoil.XData = results.coordinates.tx;
                                SW.fig_airfoil.YData = results.coordinates.tz;
                            end
                            [th, bh] = fast_thickness(results);
                            SW.fig_airfoil.Parent.Title.String = ['max bh/c = ' num2str(bh) '  and   th/c = ' num2str(th)];
                            
                            
                            if ~isempty(results.ap.bldata)
                                aoa_plot = results.ap.alpha_cl(1.4);
                                if isnan(aoa_plot)
                                    aoa_plot = results.ap.alpha_cl_max_local_limited;
                                end
                                if isempty(SW.fig_pressure)
                                    SW.fig_pressure = figure(6);                                    
                                    plot(results.ap.bldata.xun_top_range, results.ap.bldata.cpx_fun_top(aoa_plot, results.ap.bldata.xun_top_range), 'Color', [0 0.4470 0.7410]); hold on;
                                    plot(results.ap.bldata.xun_bot_range, results.ap.bldata.cpx_fun_bot(aoa_plot, results.ap.bldata.xun_bot_range), 'Color', [0 0.4470 0.7410]);
                                    grid on;
                                    xlabel('x/c');ylabel('C_p = Pressure Coefficient');
                                else
                                    SW.fig_pressure.Children.Children(1).XData = results.ap.bldata.xun_top_range;
                                    SW.fig_pressure.Children.Children(1).YData = results.ap.bldata.cpx_fun_top(aoa_plot, results.ap.bldata.xun_top_range);
                                    SW.fig_pressure.Children.Children(2).XData = results.ap.bldata.xun_bot_range;
                                    SW.fig_pressure.Children.Children(2).YData = results.ap.bldata.cpx_fun_bot(aoa_plot, results.ap.bldata.xun_bot_range);
                                end
                                SW.fig_pressure.Children.Title.String      = ['min C_p = ' num2str(max(results.ap.bldata.cpx_fun_top(aoa_plot, results.ap.bldata.xun_top_range))) ' at \alpha = ' num2str(aoa_plot) ' and Cl = ' num2str(results.ap.cl_alpha(aoa_plot))];
                            end                            
                        else
                            disp('Failed Polar: Skip Plotting.')
                        end
                    catch ME
                        disp('Error: Failed during plot, reinitialize figure handles.')
                        disp(['Error occured in: ' mfilename()])
                        disp(['Error identifier: ' ME.identifier])
                        disp( 'Error call stack: ' );
                        for is=1:length(ME.stack); disp(ME.stack(is)); end
                        SW.fig_polar   = [];
                        SW.fig_airfoil = [];
                        SW.fig_pressure= [];
                    end
                end
            end
            
        end
        
        function  results_block_cell = analyse_airfoils_block(SW, x_block_cell)
            % Analyses a block of airfoils specified by a cell array of
            % parameters in the form expected by the shape definition
            % object SD supplied to this worker object
            % This method accepts a maximum size of airfoil block equal
            % to the number of dedicated cores. For most applciations used 
            % method analyse_airfoils_list instead.
            
            % Write Airfoil Coordinate Files of each airfoil in block, in
            % the right folder
            coordinates_block_cell = SW.write_airfoils_block(x_block_cell);
            % Do the same for the simulation protocols, one for each
            % airfoil as they can be dynamized with the supplied SPD object
            SW.write_simulation_protocols_block(x_block_cell)
            
            
            ap_block_cell = SW.run_polars_block(length(x_block_cell));
            
            results_block_cell = cell(size(x_block_cell));
            
            for n_core = 1:length(x_block_cell)
                results_block_cell{n_core}.coordinates          = coordinates_block_cell{n_core};
                results_block_cell{n_core}.ap                   = ap_block_cell{n_core};
                results_block_cell{n_core}.x                    = x_block_cell{n_core};
                results_block_cell{n_core}.multisim_id_value    = SW.multisim_id_value;
                SW.output_figures(results_block_cell{n_core});
            end
            
            try 
                drawnow
            catch
                disp('Matlab is having problems plotting')
            end
        end
        
        function  results_block_cell = analyse_airfoils_block_with_polar_conditions(SW, x_block_cell, polar_conditions_block_cell)
            % Analyses a block of airfoils specified by a cell array of
            % parameters in the form expected by the shape definition
            % object SD supplied to this worker object
            % This method accepts a maximum size of airfoil block equal
            % to the number of dedicated cores. For most applciations used 
            % method analyse_airfoils_list instead.
            
            % Write Airfoil Coordinate Files of each airfoil in block, in
            % the right folder
            coordinates_block_cell = SW.write_airfoils_block(x_block_cell);
            % Do the same for the simulation protocols, one for each
            % airfoil as they can be dynamized with the supplied SPD object
            SW.write_simulation_protocols_block_with_polar_conditions(x_block_cell, polar_conditions_block_cell)
            
            
            ap_block_cell = SW.run_polars_block(length(x_block_cell));
            
            results_block_cell = cell(size(x_block_cell));
            
            for n_core = 1:length(x_block_cell)
                results_block_cell{n_core}.coordinates          = coordinates_block_cell{n_core};
                results_block_cell{n_core}.ap                   = ap_block_cell{n_core};
                results_block_cell{n_core}.x                    = x_block_cell{n_core};
                results_block_cell{n_core}.multisim_id_value    = SW.multisim_id_value;
                SW.output_figures(results_block_cell{n_core});
            end
            
            try 
                drawnow
            catch
                disp('Matlab is having problems plotting')
            end
        end
        
        function  ap_block_cell = run_polars_block(SW, varargin)
            % varargin contains number of cores to run when we don't
            % fallback to default which is all of them. This is useful when
            % we have less computations to do than cells. 
            % For example, to analyse a single airfoil, we set varargin = 1
            if ~isempty(varargin)
                varargin = varargin{1};
            end
            
            % Clean folders before execution
            SW.clean_folders(varargin)
            
%             % Write simulation protocol R/Xfoil commands
%             SW.simul_protocol.write_command_file();
            
            % Write Shell script to execute
            SW.write_shell_script(varargin);
           
            % Make computations
            SW.execute_shell_script(varargin);                              %WINDOWS LINUX MAC I think this is a correction!
            
            % Read results into polar objects (raw_data also stored in there)
            ap_block_cell = SW.load_polars(varargin);            
        end
        
        function  ap = run_polar_on_airfoil_file(SW, airfoil_file)
            % Copy Airfoil File to temporary directory for core 1 (construct a shell
            % instruction and execute it)            
            
            % Write Rfoil instruction file without dynamizer
            SW.simul_protocol.write_command_file(SW.SC.core_subdirs{1});
            
            system(['cp ' , airfoil_file  , ' ' , ...
                SW.SC.tmp_subdir, SW.SC.core_subdirs{1} , SW.SC.airfoil_filename]);
            
            % Run polar on core 1
            ap_block_cell = SW.run_polars_block(1);
            
            % Extract Aerodynamic polar object from cell array and return!
            ap = ap_block_cell{1};
        end
        
        function  results_list = analyse_airfoils_list(SW, x_list)
            % Analyses a list of airfoils specified by a cell array of
            % parameters in the form expected by the shape definition
            % object SD supplied to this worker object
            
            % Preallocate memory for cell array of results
            results_list = cell(size(x_list));
            
            % Determine into how many parallel execution blocks the list
            % should be broken down
            N_list = length(x_list);
            if SW.parallelized > 0
                c = SW.SC.N_cores;  % Use cores as defined in system context object SC
            else
                c = 1;              % Set cores to be used equal to one in case paralel execution is disabled
            end
            N_blocks = ceil(N_list/c);
            
            % Now segment execution into parallel blocks
            for n_block = 1:N_blocks
                % Determine boundaries of current execution block
                n_block_begin = (n_block - 1) * c + 1;
                n_block_end   = min(n_block_begin + c - 1, N_list);
                
                % Make cell array of parameters for this execution block
                x_block_cell = x_list(n_block_begin:n_block_end);
                
                % Execute!
                results_block_cell = SW.analyse_airfoils_block(x_block_cell);
                
                % Store Execution Results
                results_list(n_block_begin:n_block_end) = results_block_cell;                    
            end
            
            % And we are finished!
        end
        
        function  results_list = analyse_airfoils_list_with_polar_conditions(SW, x_list, polar_conditions_list)
            % Analyses a list of airfoils specified by a cell array of
            % parameters in the form expected by the shape definition
            % object SD supplied to this worker object
            
            % Preallocate memory for cell array of results
            results_list = cell(size(x_list));
            
            % Determine into how many parallel execution blocks the list
            % should be broken down
            N_list = length(x_list);
            if SW.parallelized > 0
                c = SW.SC.N_cores;  % Use cores as defined in system context object SC
            else
                c = 1;              % Set cores to be used equal to one in case paralel execution is disabled
            end
            N_blocks = ceil(N_list/c);
            
            % Now segment execution into parallel blocks
            for n_block = 1:N_blocks
                % Determine boundaries of current execution block
                n_block_begin = (n_block - 1) * c + 1;
                n_block_end   = min(n_block_begin + c - 1, N_list);
                
                % Make cell array of parameters for this execution block
                x_block_cell                = x_list(n_block_begin:n_block_end);
                polar_conditions_block_cell = polar_conditions_list(n_block_begin:n_block_end);
                
                % Execute!
                results_block_cell = SW.analyse_airfoils_block_with_polar_conditions(x_block_cell, polar_conditions_block_cell);
                
                % Store Execution Results
                results_list(n_block_begin:n_block_end) = results_block_cell;                    
            end
            % And we are finished!
        end
        
        % Philosophy for making ML foil experimental pool evaluation
        % 1. Make a generalized version of:
        %           results_list = analyse_airfoils_list(SW, x_list)
        %       ->  results_list = analyse_airfoils_list_with_polar_conditions(SW, x_list, polar_conditions_cell)
        % 2. To reach point 1, extend:
        %           results_block_cell = analyse_airfoils_block(SW, x_block_cell)
        %       ->  results_block_cell = analyse_airfoils_block_with_polar_conditions(SW, x_block_cell, polar_conditions_block_cell)
        % 3. To reach point 2, extend:
        %           write_simulation_protocols_block(SW, x_block_cell)
        %       ->  write_simulation_protocols_block_with_polar_conditions(SW, x_block_cell, polar_conditions_block_cell)
        
        
    end
end
