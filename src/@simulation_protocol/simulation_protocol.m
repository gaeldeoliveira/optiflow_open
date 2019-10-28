classdef simulation_protocol < handle
    %SIMULATION_PROTOCOL Simulation Protocol Objects define all relevant
    %procedures for simulation of an airfoil.     
    %   Fundamental Properties are listed below:
    %
    % Property: target_application              Default: 'rfoilsuc'
    % String containing calling name for target application. THis affects
    % the command file that is written, as 5.4 family codes have a slightly
    % different interface than 6.x family.
    % Valid options for command generation are:
    %       
    %       
    %
    % Property: suction_distribution_type       Default: 'none'
    % Accepts three types of suction definition:
    %   null    - no suction is applied
    %   fixed   - suction is applied based on interpolation of a
    %             suction vector
    %   'fixed_prewritten_file' -     
    %   coupled_inviscid - a material and duct underpressure is supplied, for
    %             determination of actual suction distribution at each
    %             angle of attack / speed, based on inviscid airfoil Cp(x)
    %
    %   coupled_viscous - a material and duct underpressure is supplied, for
    %             determination of actual suction distribution at each
    %             angle of attack / speed, based on viscous Cp(x)
    %
    %
    %
    % Property: geometry_correction_routine     Default: @(x) x
    % Function handle to function correcting geometrical shape parameters
    % for simulation specific shape modification. Main use is for reduction
    %  of trailing edge bluntness for greater accuracy in Xfoil. Default
    %  makes no modifications on data.    
    % Property:suction_distribution_parameters  Default: Empty struct
    % Structure containing suction distribution parameters
    %
    % 
    %
    % Property: discretization_type             Defaul: 'fixed'
    % Three possible types:
    %   fixed    - discretization does not vary with context
    %   adaptive - requests separation point determination and refines
    %              nearby
    % 
    % Property:discretization_parameters      Default: Empty struct
    % Structure containing suction discretization parameters        
    %
    % 
    % Property:operation
    % String Specifying which operation should be conducted, with
    % parameters specified in operation_parameters    
    % Possibilities are:    
    %   'alfa_polar' - calculates an angle of attack polar. 
    %                  In this case polar_properties has the following
    %                  structure:
    %                      polar_properties(1) = Lowest  Angle of Attack
    %                      polar_properties(2) = Highest Angle of Attack
    %                      polar_properties(3) = Angles of Attack increment
    %   'alfa_polar_ref_start' - calculates an angle of attack polar in two
    %    streaks starting from a reference point (with reinitialization in
    %    between streaks)
    %                      polar_properties(1) = Lowest  Angle of Attack
    %                      polar_properties(2) = Highest Angle of Attack
    %                      polar_properties(3) = Angles of Attack increment
    %                      polar_properties(4) = reference AOA
    
    %   'cl_polar'   - calculates cl polar
    %                  In this case polar_properties has the following
    %                  structure:
    %                      polar_properties(1) = Lowest  Cl
    %                      polar_properties(2) = Highest Cl
    %                      polar_properties(3) = Cl increment
    %   'cp_x'       - calculates Cp(x) at specified angles of attack. 
    %                  In this case polar_properties is a list of angles of
    %                  attack at which to determine Cp(x)
    
   
    % Command File Index
    %   1. Load Data from file    
    %   2. Repanel if needed
    %   3. Specify c/r
    %   4. Specify Suction
    %   5. Specify Re    
    %   6. Specify ITER
    %   7. Specify Output Files
    %   8. Ask for Operations (Polar or single point / BL parameters / Cpx)
    %   9. Save and Close
    
    properties
        target_application = 'RBINKOLIVEIRA_V2';
        edition = 'nsira';
        prefix = '';
        geometry_correction_routine = @(x) x;
        
        suction_distribution_type = 'none';
        suction_distribution_parameters = struct();
        
        discretization_type = 'fixed';
        discretization_parameters = struct();
        
        save_bl_parameters = 'false'
        % If true, boundary layer parameters are saved (currently only
        % under rfoilsucblind)
        
        cr = 0;
        cr_factor = 2/3;        % c/r correction factor
        
        Re = 1.6e6;             % Reynolds Number (NaN means inviscid)
        N_crit = 9;             % Critical Amplification Factor
        ITER = 300;             % Iterations
        M = 0                   % Mach number
        
        INGE = false;           % If true, toggles flag for van ingen non-linear perturbation growth
        
        N_airfoil_points = 160  % Number of airfoil points in file submited to R/Xfoil
       
        xtr_top = 1.0           % x/c of forced transition point for top side (1 = free)
        xtr_bot = 1.0           % x/c of forced transition point for bottom side (1 = free)
        
        blowing_transition = 'none'
        xtr2_top = 1.001        % x/c of transition by blowing for top side (1.001 = free)
        xtr2_bot = 1.001        % x/c of transition by blowing for bottom side (1.001 = free)
                
        repanel = 'false'       % repanel = 'false' or 'true'        
        
        operation = 'alfa_polar'
        operation_parameters = [0 20 1];
        
        EPSV = 0;               %  EPS1 convergence criterium  (only applicable on nsira edition) (EPSV = 0 keeps everything in default 1e-4 on most versions... double precision allows to go to 1e-7/1e-8 very safely and eventually to 1e-12)

        fid                     % fid pointer of currently open file
        
        SC                      % Handle to system_context object        
        airfoil_filename = 'xfin'
        SPR                     % Handle to simulation_profile object
 
    end
    
    properties(SetAccess = private)
        id
        name        
    end
        
    methods
        function simul = simulation_protocol(name , SC)
            % Constructor for simulation_protocol class object.
            simul.id = now();
            simul.fid = 1;
            
            simul.name = name;
            simul.SC   = SC;          
        end        
        function write_command_file(simul , varargin)
            % Open file for R/Xfoil commands
            % simul.fid = fopen([simul.work_dir,  simul.tmp_subdir  , 'session_tmp'],'wt');
            
            % Varargin is used to specify a subdirectory within the tmp
            % directory, for example to write the simulation protocol to a
            % core subdirectory, when in parallelized mode.
            if isempty(varargin)
                simul.fid = fopen([simul.SC.tmp_subdir , simul.name], 'wt');
                core_subdir=[];
            else
                core_subdir = varargin{1};
                simul.fid = fopen([simul.SC.tmp_subdir , core_subdir , simul.name], 'wt');
            end
            
            % Write Chapter to load airfoil coordinates
            simul.write_chapter1_load(core_subdir);
            
            % Repanel in Xfoil if activated
            if strcmp(simul.repanel ,'true')
                simul.write_chapter2_repanel();
            end
            
            % Enter OPER menu , specify Re and M numbers and max iterations
            % per point
            simul.write_chapter3_oper();
            
            % Set transition parameters (and other VPAR options if
            % applicable)
            simul.write_chapter4_bl();
            
            % Write Simulation Profile
            if ~isempty(simul.SPR)
                simul.write_chapter4_mset();
            end

            % Set suction distribution if activated
            if ~strcmp(simul.suction_distribution_type, 'none')
                simul.write_chapter5_suction();
            end
            
            % Set up Accumulation settings for later output
            simul.write_chapter6_output(core_subdir);
            
            % Perform Actual Operations
            simul.write_chapter7_oper();
            
            % Save (if not yet done) and close!
            simul.write_chapter8_save_close(core_subdir);
            
            % We're done with the command file! Close it!
            if simul.fid ~= 1
                % Only if we are writting to a file, and not to stdout
                fclose(simul.fid);
            end
        end        
        function write_chapter1_load(simul , varargin)
            if isempty(varargin)
                core_subdir=[];     
            else
                core_subdir = varargin{1};
            end
            
            fprintf(simul.fid, 'LOAD \n');
            if isempty(core_subdir)
                fprintf(simul.fid,[simul.SC.tmp_subdir  simul.SC.fs_sprinf_ap simul.airfoil_filename '\n \n']);   % New multiOS version
            else
                fprintf(simul.fid,[simul.SC.tmp_subdir simul.SC.fs_sprinf_ap]);
                fprintf(simul.fid,[core_subdir simul.SC.fs_sprinf_ap simul.airfoil_filename '\n \n']);   % New multiOS version
            end

%            fprintf(simul.fid,[ simul.airfoil_filename , ' \n \n']);       % Deprecated for greater robustness on RBINKOLIVEIRA, Fab 2014
%             if strcmp(simul.target_application(1:min(end, 5)) , 'RBINKOLIVEIRA_V2')
%                 % Load twice as first time does not always work but second does
% 
%                 
%                 fprintf(simul.fid, 'LOAD \n');
%                 fprintf(simul.fid, [ simul.airfoil_filename , ' \n \n']);
%             else
%                 %Load Profile Points
%                 fprintf(simul.fid, ['LOAD ' , [simul.SC.tmp_subdir core_subdir simul.airfoil_filename], ' \n \n']);
%                 fprintf(simul.fid, 'PLOP \n');    % Set Xplot window
%                 fprintf(simul.fid, 'W 0.4\n');    % Reduce size to 40% fo screen
%                 fprintf(simul.fid, 'G \n');       % Disable, if version allows
%                 fprintf(simul.fid, '\n');         % Same as above
%             end

        end
        function write_chapter2_repanel(simul)
            if strcmp(simul.reinterpolate, 'true')
                % Spline reinterpolation is only applicable to xfoil 6.x  and
                % hence does not work on rfoil(suc) based on xfoil 5.4
                if strcmp(simul.target_application(1:min(end, 5)) , 'xfoil')
                    fprintf(simul.fid,'GDES\n');          % Enter Data treatment menu
                    fprintf(simul.fid,'CADD\n');          % Interpolate Points with Spline For smooth Leading Edge
                    fprintf(simul.fid,'\n \n \n \n \n');  % Let X foil find good values
                else
                    disp(['Warning: Attempted Reinterpolation on incompatible target application (' , simul.target_application, ')'])
                    % Abd al Malik - C'Est Du Lourd!
                end
            end            
            % Emily Jane White - Victorian America            
            fprintf(simul.fid,'PANE\n');          % Make Panelling
        end
        function write_chapter3_oper(simul)
            fprintf(simul.fid, 'OPER\n');         % Enter oper menu            
            % If viscosity on activate and specify Reynolds Number
            if ~isnan(simul.Re)
                fprintf(simul.fid, 'VISC \n');
                fprintf(simul.fid, '%0.4e \n' , simul.Re); 
            end
            
            if strcmp(simul.target_application(1:min(end, 5)) , 'xfoil')
                fprintf(simul.fid, 'M \n %0.4e \n'    , simul.M);     % Set Mach Number
                fprintf(simul.fid, 'ITER \n %0.4i \n' , simul.ITER);  % Set Maximum Iteration limit per point
            else
                fprintf(simul.fid, 'MACH \n %0.4e \n' , simul.M);     % Set Mach Number
            end            
        end
        function write_chapter4_bl(simul)
            % Enter VPAR menu
            fprintf(simul.fid, 'VPAR\n');
            
            % Set Critical Amplification Factor            
            fprintf(simul.fid, 'N \n %0.4e \n', simul.N_crit);
            
            % Toggle van Ingen non-linear growth flag, if applicable
            if (simul.INGE == true)
                fprintf(simul.fid, 'INGE \n');
            end
            
            % Set Transition Points (xtr = 1 means free transition)
            fprintf(simul.fid, 'XTR \n %0.4e \n %0.4e \n', simul.xtr_top , simul.xtr_bot);
            
            % % Deprecated: Change was embedded back into into Rfoil
            % % BY RICARDO !! Change flags to match RFOIL of Nando
            % fprintf(simul.fid, 'l_sc \n');
            % fprintf(simul.fid, 'l_cf \n');
            
            % command for transition by blowing (only in rfoil)
            if ~strcmp(simul.target_application(1:min(end, 5)) , 'xfoil') && strcmp(simul.blowing_transition , 'on')
                fprintf(simul.fid,'XTR2 \n %0.4e \n %0.4e \n', simul.xtr2_top , simul.xtr2_bot);
            end
               
             if strcmp(simul.edition , 'nsira')
                 if simul.EPSV ~= 0
                     fprintf(simul.fid, 'EPSV \n');
                     fprintf(simul.fid, [num2str(simul.EPSV) , ' \n']);
                 end
             end
            
            % Come back to OPER menu
            fprintf(simul.fid, '\n');
        end
        function write_chapter4_mset(simul)
            % Enter VPAR menu
            fprintf(simul.fid, 'MSET\n');
            
            % Write Bernstein Polynomial Coefficients of Cf Parametrization
            fprintf(simul.fid, 'TA1C \n %0.6e \n', simul.SPR.TA1C );
            fprintf(simul.fid, 'TA2C \n %0.6e \n', simul.SPR.TA2C );
            fprintf(simul.fid, 'TA3C \n %0.6e \n', simul.SPR.TA3C );
            fprintf(simul.fid, 'TA4C \n %0.6e \n', simul.SPR.TA4C );
            fprintf(simul.fid, 'TA5C \n %0.6e \n', simul.SPR.TA5C );
            fprintf(simul.fid, 'TA6C \n %0.6e \n', simul.SPR.TA6C );
            
            % Write Bernstein Polynomial Coefficients of Cf Parametrization
            fprintf(simul.fid, 'TA1H \n %0.6e \n', simul.SPR.TA1H );
            fprintf(simul.fid, 'TA2H \n %0.6e \n', simul.SPR.TA2H );
            fprintf(simul.fid, 'TA3H \n %0.6e \n', simul.SPR.TA3H );
            fprintf(simul.fid, 'TA4H \n %0.6e \n', simul.SPR.TA4H );
            fprintf(simul.fid, 'TA5H \n %0.6e \n', simul.SPR.TA5H );
            fprintf(simul.fid, 'TA6H \n %0.6e \n', simul.SPR.TA6H );
            
            % Write Bernstein Polynomial Coefficients of Cf Parametrization
            fprintf(simul.fid, 'THMN \n %0.6e \n', simul.SPR.THMIN);
            fprintf(simul.fid, 'THMX \n %0.6e \n', simul.SPR.THMAX);
            fprintf(simul.fid, 'TDCF \n %0.6e \n', simul.SPR.TDCF );
            
            fprintf(simul.fid, 'DIAG \n');
            
            % Come back to OPER menu
            fprintf(simul.fid, '\n');
        end
        function write_chapter5_suction(simul)
            if strcmp(simul.suction_distribution_type,'fixed_prewritten_file')
                % Enter Suction Design Menu
                fprintf(simul.fid, 'VDES\n');
                % Load Suction Distribution from file
                fprintf(simul.fid, ['LOAD\n' ,  simul.SC.suc_filename , '\n']);
                % Set as active and return to OPER menu
                fprintf(simul.fid, 'EXEC\n\n');
            end
            
            if strcmp(simul.suction_distribution_type,'target_H')
                
                % Enter Suction Design Menu
                fprintf(simul.fid, 'VDES\n');
                
                % Specify Suction Distribution Parameters
                fprintf(simul.fid, 'SPAR\n');
                fprintf(simul.fid, 'TMET 1\n');   % TMET=1 -> Generate for target H
                fprintf(simul.fid, 'CHAL \n');    % Change all points at the same time (more robust)                
                
                htur = simul.suction_distribution_parameters.htur;
                fprintf(simul.fid, ['HTUR ', num2str(htur) , '\n']); % Set target suction distribution
                
                fprintf(simul.fid, '\n');         % Return to VDES menu
                
                % Mark Suction Area (should be customizable later)
                mark_first = num2str(simul.suction_distribution_parameters.mark(1));
                mark_last  = num2str(simul.suction_distribution_parameters.mark(2));
                fprintf(simul.fid, ['MARK ' , mark_first , ' ' , mark_last , ' \n']);
                
                % Return to OPER menu
                fprintf(simul.fid, '\n');
                
            end
            
            if strcmp(simul.suction_distribution_type,'rfoil_NSM2')
                % Enter Suction Design Menu
                fprintf(simul.fid, 'VDES\n');
                
                % Specify Suction Mode to 2 (parametric predefined suction)
                fprintf(simul.fid, 'NSM \n');
                fprintf(simul.fid, '2\n');
                
                % Extract values from suction parameters structure
                ssuc = simul.suction_distribution_parameters.ssuc;
                xsuc = simul.suction_distribution_parameters.xsuc;
                vsuc = simul.suction_distribution_parameters.vsuc;
                
                % Set suction side and interval
                fprintf(simul.fid, 'ssuc \n');
                fprintf(simul.fid, [num2str(ssuc) , '\n']);
                fprintf(simul.fid, 'xsuc \n');
                fprintf(simul.fid, [num2str(xsuc(1)) , '\n']);
                fprintf(simul.fid, [num2str(xsuc(2)) , '\n']);
                
                % Set suction speed
                
                fprintf(simul.fid, 'vsuc \n');
                fprintf(simul.fid, [num2str(vsuc) , '\n']);
                
                % Return to OPER menu
                fprintf(simul.fid, '\n');
            end
            
        end
        function write_chapter6_output(simul , varargin)
            if isempty(varargin)
                core_subdir=[];
            else
                core_subdir = varargin{1};
            end
            
            
            if ~strcmp(simul.target_application(1:min(end, 5)) , 'xfoil')
                % If Boundary Layer Parameters are to be saved activate
                % Auto-fill
                if strcmp(simul.save_bl_parameters, 'true')
                    fprintf(simul.fid, 'ARBL \n');
                    % NOTE: This must be done before PACC is called, as
                    % otherwise the saved file is very, very strange
                end
            end            
           
            % Activate Polar accumulation mode
            fprintf(simul.fid, 'PACC\n');
            % Specify Polar accumulation file
            fprintf(simul.fid, [simul.SC.tmp_subdir simul.SC.fs_sprinf_ap core_subdir simul.SC.fs_sprinf_ap simul.SC.polar_filename , '\n']);
            % Specify Polar dump file
            fprintf(simul.fid, [simul.SC.tmp_subdir simul.SC.fs_sprinf_ap core_subdir simul.SC.fs_sprinf_ap simul.SC.polar_dump_filename , '\n']);
            
            % Mute convergence output
            fprintf(simul.fid, 'MUTE\n');
            
        end    
        function write_chapter7_oper(simul)
            if strcmp(simul.operation , 'alfa_polar')
                fprintf(simul.fid, 'ASEQ \n');
                fprintf(simul.fid, [num2str(simul.operation_parameters(1)) , ' \n']);
                fprintf(simul.fid, [num2str(simul.operation_parameters(2)) , ' \n']);
                fprintf(simul.fid, [num2str(simul.operation_parameters(3)) , ' \n']);
            end
            
            
            % If applicable run some points before storing (used to prevent
            % convergence towards undesirable solution
            if strcmp(simul.operation , 'alfa_polar_ref_start')
                
                % Compute first streak, from ref_start to end
                
                fprintf(simul.fid, 'ASEQ \n');
                fprintf(simul.fid, [num2str(simul.operation_parameters(4)) , ' \n']);
                fprintf(simul.fid, [num2str(simul.operation_parameters(2)) , ' \n']);
                fprintf(simul.fid, [num2str(simul.operation_parameters(3)) , ' \n']);
                
                % Reinitialize boundary layer
                fprintf(simul.fid, 'VPAR \n');          % Go into VPAR menu for Rfoil compatibility
                fprintf(simul.fid, 'INIT \n');          % Reinitialize
                fprintf(simul.fid, '\n');               % Return to OPER menu
                
                % Compute second streak, from ref_start to beginning                
                fprintf(simul.fid, 'ASEQ \n');
                fprintf(simul.fid, [num2str(simul.operation_parameters(4)- simul.operation_parameters(3)) , ' \n']);
                fprintf(simul.fid, [num2str(simul.operation_parameters(1)) , ' \n']);
                fprintf(simul.fid, [num2str( - simul.operation_parameters(3)) , ' \n']);
            end
            
            
            if strcmp(simul.operation , 'cl_polar')
                fprintf(simul.fid, 'CSEQ \n');
                fprintf(simul.fid, [num2str(simul.operation_parameters(1)) , ' \n']);
                fprintf(simul.fid, [num2str(simul.operation_parameters(2)) , ' \n']);
                fprintf(simul.fid, [num2str(simul.operation_parameters(3)) , ' \n']);
            end            
            if strcmp(simul.operation , 'cp_x')
                for n_alfa = 1:length(simul.operation_parameters)
                    fprintf(simul.fid, ['alfa ' , num2str(simul.operation_parameters(n_alfa)) , ' \n']);
                    fprintf(simul.fid, ['CPWR ' , simul.SC.svcp_filename , num2str(alfa) , ' \n']);                    
                end                                                    
            end
                        
            if strcmp(simul.operation , 'target_H_alfa_polar')
                % The polar must be scripted point by point in target_H
                % mode.
                % If PACC is enabled, the points calculated in the suction 
                % design routine are accumulated and storage arrays will 
                % overflow very fast... So it is better to calculate
                % target_H with polar accumulation disabled and then store
                % results
                
                % Make list of AOA's at which to calculate suction
                % distribution
                alfa_vector = simul.operation_parameters(1):simul.operation_parameters(3):simul.operation_parameters(2);
                N_alpha = length(alfa_vector);
  
                for n_alpha = 1:N_alpha                                                            
                    % For AOA's below maximum AOA for suction generation 
                    if alfa_vector(n_alpha) < simul.suction_distribution_parameters.AOAmax
                        fprintf(simul.fid, 'ITER \n %0.4i \n' , 10);  % Set Maximum Iteration limit to 10 to speed up suction computation
                        
                        %Generate Suction Distribution
                        fprintf(simul.fid, 'PACC \n');      % Disable PACC
                        
                        % First Generate Initial Point Without suction for better
                        % base estimation
                        if n_alpha == 1
                            fprintf(simul.fid, ['ALFA ' , num2str(alfa_vector(n_alpha)) , '\n']);
                        end                        
                        
                        fprintf(simul.fid, 'VDES \n');      % Enter VDES menu
                        
                        
                        % Generate Scaled Base Initial Suction for first point, 
                        % for faster convergence better
                        if n_alpha == 1
                            fprintf(simul.fid, 'BASE\n');
                            fprintf(simul.fid, 'SCAL 0.3\n');
                            fprintf(simul.fid, 'EXEC\n');
                        end
                        
                        
                        fprintf(simul.fid, ['ALFA ' , num2str(alfa_vector(n_alpha)) , ' \n']);      % Actually generate Suction Distribution for this point
                        
                        % Allow disabling for Polar storage for
                        % comparison/validation purposes
                        if simul.suction_distribution_parameters.disabled == 1
                            fprintf(simul.fid, 'SCAL 0\n');
                            fprintf(simul.fid, 'EXEC 0\n');
                        end
                                                
                        fprintf(simul.fid, '\n');           % Return to OPER menu
                        
                        % Re-Activate Polar accumulation mode
                        fprintf(simul.fid, 'PACC\n');
                        fprintf(simul.fid, [ simul.SC.polar_filename , '\n']);     % Specify Polar accumulation file
                        fprintf(simul.fid, [ simul.SC.polar_dump_filename , '\n']);    % Specify Polar dump file                        
                    end
                    
                    % Compute Point for accumulation
                    fprintf(simul.fid, 'ITER \n %0.4i \n' , simul.ITER);  % Set Maximum Iteration limit for accumulation
                    fprintf(simul.fid, ['ALFA ' , num2str(alfa_vector(n_alpha)) , ' \n']);                    
                end 
            end            
        end
        function write_chapter8_save_close(simul, varargin)
            % WRITE_CHAPTER_9_SAVE_CLOSE
            % Writes final part of R/Xfoil commandfile according to options            
            if isempty(varargin)
                core_subdir=[];
            else
                core_subdir = varargin{1};
            end
            
            % If Boundary Layer Parameters are to be saved activate
            % Auto-fill (only under rfoilsucblind)
            if ~strcmp(simul.target_application(1:min(end, 5)) , 'xfoil')
                if  strcmp(simul.save_bl_parameters, 'true')
                    % For this:
                    %   not to crash everything we need less than 150-200 alpha
                    %   to work we need very few points (explore this better
                    
                    
                    % Request saving of BL parameters
                    fprintf(simul.fid, 'SVBL \n');
                    % At all stances (0 would produce a SILANT file)
                    fprintf(simul.fid, '1 \n');
                    % Specify file
                    fprintf(simul.fid , [simul.SC.tmp_subdir simul.SC.fs_sprinf_ap core_subdir simul.SC.fs_sprinf_ap simul.SC.svbl_filename , '\n']);
                    %fprintf(simul.fid, [simul.SC.tmp_subdir simul.SC.fs_sprinf_ap core_subdir simul.SC.fs_sprinf_ap simul.SC.polar_filename , '\n']);
                end            
            end
            
            fprintf(simul.fid, '\n\n');                     % Return to root menu
            fprintf(simul.fid, 'QUIT\n');                   % Finish execution
                   
        end
        
       
    end
    
end

