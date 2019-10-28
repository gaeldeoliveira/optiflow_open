classdef swafford_profile < handle
    %SWAFFORD_PROFILE is a class to facilitate the evaluation of swafford
    %   profiles while minimizing computational overhead in the particular
    %   case of integrations along the y coordinate.
    %   Meant to help development of a closure relation for the (Plasma
    %   Force Field) Energy Interaction Coefficient (see Report!)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Integral Boundary Layer Integrator (ODE)
    %           Plasma Development Tool
    %
    %   August 2014, GNU-GPLv3 or later
    %   Gael de Oliveira, Ricardo Pereira
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Swafford velocity profile generator
    %          following:
    %              Analytical Approximation of Two-Dimensional Separated
    %              Turbulent Boundary-Layer Velocity Profiles
    %              Swafford, T.W., AIAA Journal Vol.21, N.6, June 1983
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties
        % Primary Inputs
        hk  = 0;    % Kinematic Shape Factor (delta* / theta)
        rt  = 0;    % Reynolds Theta (Momentum Thickness Reynolds Number)
        msq = 0;    % Mach Number (for data consistency in Cf calls only!)
        
        % Dependent Variables (on h and rt)
        cf          % Step 2
        s           % Step 2
        ue_p        % Step 3
        u_ue_2      % Step 4
        u_ue_5      % Step 5
        g_2         % Step 6
        g_5         % Step 7
        b           % Step 8
        a           % Step 9
        
        % 
        
        
    end
    
    methods
        function SP = swafford_profile()
            % Constructor Method
        end
        
        function update_hk_rt_pair(SP, hk, rt)
            % Check wether new h and rt pair is different from original
            if or(or(not(SP.hk==hk) , not(SP.rt==rt)),isempty(SP.cf))
                % If they are different, update
                SP.rt = rt;
                SP.hk = hk;
                
                % And also update dependencies
                SP.update_dependent_variables();                
            end
            % Otherwise, do nothing!
        end
        
        function update_dependent_variables(SP)
                        
            % % Step 2
            % Calculate Skin Friction
            % Static method invocation!
            %SP.cf= cft( SP.hk, SP.rt, SP.msq);                     % Original Rfoil Cf relation, (innapropriate for profile generation)
            SP.cf = swafford_profile.cft_swafford(SP.hk, SP.rt);    % Swafford profile generation Cf closure. See method description for more details on why this is the right choice!
                        
            % Calculate ue_p = ue+
            SP.ue_p = sqrt(2 ./ abs(SP.cf));
            
            % % Step 3
            % Compute Sign of Skin Friction
            % s = cf / abs(cf)
            SP.s = sign(SP.cf);
            
            % % Step 4
            % Compute u over ue (2)
            SP.u_ue_2 = (1 ./ 1.95) * ( atanh((8.5-SP.hk)/7.5) - 0.364 );
            
            % % Step 5
            % Compute u over ue (5)
            SP.u_ue_5 = 0.155 + 0.795 * sech(0.51 * (SP.hk-1.95));
            
            % % Step 6
            % Compute g(2)
            SP.g_2 = (SP.u_ue_2 - (SP.s ./ (0.09 * SP.ue_p)) .* atan(0.18 * SP.rt./SP.ue_p)) ...
                ./ (1 - (SP.s * pi()) ./ (0.18*SP.ue_p) );
            
            % % Step 7
            % Compute g(5)
            SP.g_5 = (SP.u_ue_5 - (SP.s ./ (0.09 * SP.ue_p)) .* atan(0.45 * SP.rt./SP.ue_p)) ...
                ./ (1 - (SP.s * pi()) ./ (0.18*SP.ue_p) );
            
            % % Step 8
            % Compute b (recall that in matlab log == ln)
            SP.b = log( atanh((SP.g_2).^2) ./ atanh((SP.g_5).^2) ) ./ log(2/5);
            
            % % Step 9
            % Compute a
            SP.a = atanh(SP.g_2.^2) / (2.^(SP.b));
            
        end
        
        function [u_over_ue, u_p]= evaluate_profile(SP, y_over_theta)
            % % Step 10
            % Evaluate profile!
            y_p = (SP.rt ./ SP.ue_p) .* y_over_theta;
            
            u_p = (SP.s ./ 0.09) * atan(0.09 * y_p) + (SP.ue_p - SP.s*pi()/0.18) .* sqrt( tanh(SP.a .*  y_over_theta.^(SP.b)) );
            
            u_over_ue = u_p / SP.ue_p;
            
        end
        
        function u_over_ue = u_over_ue(SP, hk, rt, y_over_theta)
            % Update hk and rt dependencies
            SP.update_hk_rt_pair(hk, rt);
            
            % Evaluate Profile
            u_over_ue = SP.evaluate_profile(y_over_theta);
        end
        
        function plot(SP, hk, rt, varargin)
            % Update hk and rt dependencies
            SP.update_hk_rt_pair(hk, rt);
            
            y_over_theta_vector = linspace(0,10,100);
            u_over_ue_vector = SP.evaluate_profile(y_over_theta_vector);
            plot(u_over_ue_vector, y_over_theta_vector, varargin{:})
            grid on
            xlabel('U over Ue')
            ylabel('y over theta')   
        end
        
    end
    
    methods(Static)
        function cf = cft_swafford(hk, rt)
            % Swafford Cf function, as written in the paper, needed to
            % generate profiles with the expected shape factor (H). 
            %
            % If another Cf function is used, like Rfoil's cf function,
            % another shape factor (H) appears upon integration, thereby
            % creating inconsistencies in the derivation of closure
            % relations. 
            %
            % For example, for [hk = 1.8 ; rt = 2000] and using rfoil's CFT
            % we recover a hk = 1.7776 upon integration wi 
            
            % Break up computation in two parcels
            first_parcel  = 0.3 * exp(-1.33 * hk) ./ (log10(rt).^(1.74+0.31*hk));
            second_parcel = (1.1e-4) * (tanh(4 - hk./0.875) - 1);
            
            % Sum to obtain Swafford's Cf relation, from Whitfield (I think!)
            cf = first_parcel + second_parcel;            
        end
    end
    
end