classdef simulation_protocol_dynamizer < handle
    %SIMULATION_PROTOCOL_DYNAMIZER This class takes the dummy_arguments
    %from the parameter vector and uses them to modify the simulation
    %parameters, for example controlling suction.
    % Currently this class only acts on suction, with two possible modes,
    % fixed suction speed, or fixed Cq, and the optimized parameters being
    % the suction area begin and end x/c coordinates.
    % However, it is called simulation_parameter_dynamizer as the
    % framework let's this class do anything you like on the simulation!   
    %
    %
    % A creative application of this class:
    % The overall framework provides the possibility to use this class for 
    % tunning of R/Xfoil boundary parameters to sets of experimental
    % results. Obviously this requires modification of this class and a 
    % new cost function to be constructed
    
    
    properties
       suction_distribution_type = 'rfoil_NSM2'
       mode = 'Fixed_Suction_Speed' ;   
       %  Mode controls use of dummy parameters, without affecting blade
       %  localization. If can be either :
       %         'Inactive'
       %         'Fixed_Suction_Speed'
       %         'Fixed_Suction_CQ'
       
       
       fixed_Cqm  = -0.003*0.2          % Cq value to use in 'Fixed_Suction_CQ' mode. Cq is defined here as Cq = (Vs / Uinf) * (length_suction_zone / chord)       
       fixed_Cqm_vsuc_max   = -0.02     % Maximum allowable suction              
       
       fixed_vsuc = -0.003              % Suction speed value to use in 'Fixed_Suction_Speed' mode
       ssuc = 1                         % Side on which to apply suction in 'Fixed_Suction_Speed' and 'Fixed_Suction_Mass' modes  (side 1 = upper, side 2 = lower)
       
       simul_protocol                   %  Simulation Protocol affected by this class
       N_simulation_parameters          %  Number of simulation arguments       
       
       localize_profiles  = false       % Use a blade object to localize profiles based on their thickness
       
       % Optional  Object Handles nessary for profile localization       
       BD                               % Handle to blade object, if applicable
       SD                               % Shape Definition Object
       
    end
    
    methods
        function SPD = simulation_protocol_dynamizer(mode , simul_protocol , varargin)
            % Assign compulsory fields
            SPD.mode = mode;
            SPD.simul_protocol = simul_protocol;
            % Assign blade object handle if applicable
%             if ~isempty(varargin)
             if length(varargin) == 2
                SPD.BD = varargin{1};
                SPD.SD = varargin{2};
                   % When Profile Localization handles are provided, enable
                   % it by default
                SPD.localize_profiles = true;
            end
            
            % Set number of parameters assigned to simulation control
            if strcmp(mode, 'Inactive')
               SPD.N_simulation_parameters = 0;
            end
            if strcmp(mode, 'Fixed_Suction_Speed')
               SPD.N_simulation_parameters = 2; 
            end
            if strcmp(mode, 'Fixed_Suction_Mass')
               SPD.N_simulation_parameters = 2;                 
            end            
        end                
                
        function localize_profile_on_blade(SPD, x_parameters)
            % Function that takes cares of localizing the profile on the
            % blade, that is:
            %  Determine radius of profile from blade thickness distribution
            %  Set Re and cr according to blade properties
            
            
            % Find Thickness of profile (upgrade this snippet to cleaner
            % code one day!)
            [~ , tz] = SPD.SD.generate_coordinates(...
                    300 , x_parameters);
                
            tc_profile =  max(tz) - min(tz);                       
            
            % Request conditions of profile on blade from blade object            
            % Stat by local solidity
            cr = SPD.BD.cr_from_tc(tc_profile);
            % and finish with Reynolds
            Re = SPD.BD.Re_from_tc(tc_profile);
            
            % Now set modifications in simulation_protocol !
            SPD.simul_protocol.cr = cr;
            SPD.simul_protocol.Re = Re;
        end        
                        
        function modify_simulation_protocol(SPD, x_parameters)
            % Function called by other objects to introduce applicable
            % modifications to simulation protocol
            
            % First Localize Blade (if applicable)
            if SPD.localize_profiles == true 
                SPD.localize_profile_on_blade(x_parameters)
            end
            
            
            % Extract simulation parameters from overall parameters vector
            simul_parameters = x_parameters((end-SPD.N_simulation_parameters+1):end);
            
            if strcmp(SPD.mode, 'Inactive')
               % Do nothing 
            end
            if strcmp(SPD.mode, 'Fixed_Suction_Speed')
                SPD.modify_simulation_protocol_FSV(simul_parameters)
            end
            if strcmp(SPD.mode, 'Fixed_Suction_Mass')
                SPD.modify_simulation_protocol_FSQ(simul_parameters)                
            end
        end
        
        function modify_simulation_protocol_FSV(SPD , simul_parameters)
            % Build Suction distribution parameters structure according to
            % parameters!
            
            % Extract suctin parameters and apply variable change to obtain
            % beggining and end of suction region
            l_s  = simul_parameters(1);
            x2s = simul_parameters(2);
            x1s = x2s - l_s;
                        
            SPD.simul_protocol.suction_distribution_parameters.ssuc = SPD.ssuc;
            SPD.simul_protocol.suction_distribution_parameters.xsuc = [x1s x2s];
            SPD.simul_protocol.suction_distribution_parameters.vsuc = SPD.fixed_vsuc;
            disp(['Suction :    Length = ' , num2str(l_s)  , '   from : ' , num2str(x1s) , ' to ' , num2str(x2s) , '    at  vs =' , num2str(vsuc)]);
        end
        
        function modify_simulation_protocol_FSQ(SPD , simul_parameters)
            % Build Suction distribution parameters structure according to
            % parameters!
            
            % Extract suction parameters and apply variable change to obtain
            % beggining and end of suction region
            l_s  = simul_parameters(1);
            x2s = simul_parameters(2);
            x1s = x2s - l_s;
            
            SPD.simul_protocol.suction_distribution_parameters.ssuc = SPD.ssuc;
            SPD.simul_protocol.suction_distribution_parameters.xsuc = [x1s x2s];
            
            % Vs/Uinf speed now needs to be computed from Cq
            % Cq definition varies throughout literature. Here, Cq is
            % defined as Cq = (Vs / Uinf) * (length_suction_zone / chord)
            
            % Find length of suction zone (abs used to avoid 
            % Now find vsuc from suction zone length ls and Cq
            vsuc = SPD.fixed_Cqm / l_s;
            
            % Limit vsuc to reasonable values!
            vsuc = max(SPD.fixed_Cqm_vsuc_max , vsuc);
            
            % And set it into the suction parameters structure
            SPD.simul_protocol.suction_distribution_parameters.vsuc = vsuc;
            disp(['Suction   Length = ' , num2str(l_s)  , ' from : ' , num2str(x1s) , ' to ' , num2str(x2s) , '  at  vs =' , num2str(vsuc)]);
        end
    end
    
end




